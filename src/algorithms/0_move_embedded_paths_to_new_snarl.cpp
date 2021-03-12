#include "0_oo_normalize_snarls.hpp"

namespace vg {
namespace algorithms {
/*
* Goal: break move_path_to_snarl into easier-to-read, composite parts.
* 
* Make it easier to take into account "backwards" bool when wanting to iterate through
* snarl backwards.
* 
* Temp notes to self:
* for compare syntax, see https://www.geeksforgeeks.org/stdstringcompare-in-c/
* syntax 2 was what I wanted when I looked it up.
* Example Error message:     
*    if (get_path_handle_of_step(segment_begin) != get_path_handle_of_step(segment_end)) {
*        cerr << "error:[VG] attempted to rewrite segment delimited by steps on two separate paths" << endl;
*        exit(1);
*    } 
*
*/

vector<pair<vector<handle_t>, int>> SnarlNormalizer::find_possible_path_starts(const handle_t &leftmost_handle, const handle_t &rightmost_handle, const pair<bool, bool> &path_spans_left_right)
{
    vector<pair<vector<handle_t>, int>> possible_paths;

    // if path starts at leftmost handle, then path's leftmost extension is at the
    // beginning of the leftmost handle.
    if (path_spans_left_right.first)
    {
        vector<handle_t> path;
        path.push_back(leftmost_handle);
        pair<vector<handle_t>, int> path_loc = make_pair(path, _graph.get_sequence(leftmost_handle).size());
        possible_paths.push_back(path_loc);
    }
    return possible_paths;
}
/**
 * extend_possible_paths
 * 
 * @param  {vector<pair<vector<handle_t>} undefined : 
 * @param  {int>>} possible_path_starts             : 
 * @param  {string} path_str                        : 
 * @param  {handle_t} leftmost_handle               : 
 * @param  {handle_t} rightmost_handle              : 
 * @param  {pair<bool} undefined                    : 
 * @param  {bool>} path_spans_left_right            : 
 * @return {vector<handle_t>}                       : empty if no path; otherwise, the sequence of handles representing a path which matches the path_str.
 */
vector<handle_t> SnarlNormalizer::extend_possible_paths(vector<pair<vector<handle_t>, int>> &possible_path_starts, const string &path_str, const handle_t &leftmost_handle, const handle_t &rightmost_handle, const pair<bool, bool> &path_spans_left_right)
{
    // cerr << "path string (note: should be left-to-right at this point, e.g. TTACT, not AGTAA: " << path_str << endl;
    // cerr << "leftmost handle id and seq: " << _graph.get_id(leftmost_handle) << " " << _graph.get_sequence(leftmost_handle) << endl;
    vector<handle_t> correct_path;
    // Now that we have all the possible leftmost starting positions for the path in
    // possible_paths, search to the right of those positions for possible extensions
    // of the path.
    //
    // If there are two possible extensions from a single position, make two entries in
    // possible_paths.
    //
    // Continue until extending each possible path location, until we find a path that
    // reaches (and properly includes) the sink. If there is no such path, then throw an
    // exception saying that we couldn't find the path.
    int times=0;
    while (!possible_path_starts.empty() && correct_path.empty())
    {
        // cerr << "walked through while loop " << times << " times." << endl;

        times += 1;
        // take a path off of possible_path_starts, which will be copied for every iteration
        // through _graph.follow_edges, below:
        pair<vector<handle_t>, int> cur_possible_path = possible_path_starts.back();
        possible_path_starts.pop_back();

        string possible_path_str;
        for (handle_t handle : cur_possible_path.first)
        {
            // cerr << "cur_possible_step id: " << _graph.get_id(handle) << endl;
            possible_path_str += " " + _graph.get_sequence(handle);
        }

        // cerr << "cur_possible_path is size " << cur_possible_path.first.size() << " with seq length " << cur_possible_path.second << " and sequence " << possible_path_str << endl;
        // extend the path through all right-extending edges to see if any subsequent
        // paths still satisfy the requirements for being a possible_path:
        _graph.follow_edges(
            get<0>(cur_possible_path).back(), false, [&](const handle_t &next) {
                // make a copy to be extended for through each possible next handle in
                // follow edges.
                pair<vector<handle_t>, int> possible_path = cur_possible_path;

                // decide if "next" is a valid extension of possible_path.
                string next_seq = _graph.get_sequence(next);
                // cerr << "next_seq is from candidate node for extension. next_seq: " << next_seq << endl;
                //todo: test that this compare functions properly.
                if ((path_str.compare(possible_path.second, next_seq.size(), next_seq) == 0))
                {
                    // cerr << "candidate for extension passed!" << endl;
                    possible_path.first.push_back(next);
                    pair<vector<handle_t>, int> new_possible_path = make_pair(possible_path.first, possible_path.second + next_seq.size());
                    // If we've reached the end of the path, we've either reached the proper end point (anywhere if path_spans_left_right.second==false; else, the rightmost node.), or we've yet to find a valid new_possible_path.
                    if ((!path_spans_left_right.second && new_possible_path.second>= path_str.size()) 
                        || (path_spans_left_right.second && _graph.get_id(next) == _graph.get_id(rightmost_handle) && new_possible_path.second == path_str.size()))
                    {
                        // we found the correct path.
                        correct_path = new_possible_path.first;
                        return;
                    }
                    else 
                    {
                        possible_path_starts.push_back(new_possible_path);
                    }
                }
                else 
                {
                    // cerr << "candidate for extension failed." << endl;
                }
                // return false;
            });
    }
    if (correct_path.size() == 0)
    {
        cerr << "************UNIT_TEST for extend_possible_paths************" << endl;
        cerr << "no correct path in snarl found. output path size is zero." << endl;
        cerr << "************END-UNIT_TEST for extend_possible_paths.************"<< endl;
    }

    return correct_path;
}
/**
 * SnarlNormalizer::move_path_to_new_snarl 
 * 
 * @param  {pair<step_handle_t, step_handle_t>} old_path : 
 * @param  {id_t} source                                 : the source for the new snarl.
 * @param  {id_t} sink                                   : the sink for the new snarl.
 * @param  {pair<bool} undefined                         : 
 * @param  {bool>} path_spans_left_right                 : 
 * @param  {bool} path_directed_left_to_right            : 
 * @return {pair<step_handle_t, step_handle_t>}          :
 */
pair<step_handle_t, step_handle_t> SnarlNormalizer::move_path_to_new_snarl(const pair<step_handle_t, step_handle_t> & old_path, const id_t &leftmost_id, const id_t &rightmost_id, const pair<bool, bool> &path_spans_left_right, const bool &path_directed_left_to_right)
{
    /*
* This should return the series of handles, from left to right if path_left_to_right==true (else vice-versa), that the path should move to.
* 
* Or returns None if the proposed "valid_starting_index" didn't pan out to give a good
* path in following handles.
*/
    // if path doesn't span both source and sink, I need to address that. But for the
    // first iteration of this algorithm, I'll dodge that question.
    // todo: address paths that don't span both source and sink.
    if (!(path_spans_left_right.first and path_spans_left_right.second))
    {
        cerr << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
        cerr << "PATH DOESN'T SPAN SOURCE AND SINK! THIS IS CURRENTLY UNSUPPORTED. SNARL WILL BE NORMALIZED, BUT PATH WON'T BE INCLUDED." << endl;
        cerr << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
        pair<step_handle_t, step_handle_t> no_path;
        return no_path;
    }

    // get the path_string from the handles in the old_path:
    string path_str;
    step_handle_t cur_step = old_path.first;
    string path_name = _graph.get_path_name(_graph.get_path_handle_of_step(old_path.first)); // used for unit tests at bottom.
    vector<handle_t> old_path_location;
    while (cur_step != old_path.second)
    {
        // cerr << "orientation of cur_step: " << _graph.apply_orientation
        path_str += _graph.get_sequence(_graph.get_handle_of_step(cur_step));
        //todo: note following line of for loop is for debug purposes. delete?
        old_path_location.push_back(_graph.get_handle_of_step(cur_step));
        cur_step = _graph.get_next_step(cur_step);
    }

    // cerr << "path string as originally extracted: " << path_str << endl;
    handle_t leftmost_handle = _graph.get_handle(leftmost_id);
    handle_t rightmost_handle = _graph.get_handle(rightmost_id);

    // make path_str read left-to-right.
    if (!path_directed_left_to_right)
    {
        path_str = reverse_complement(path_str);
    }
        
    //     leftmost_handle = _graph.get_handle(sink);
    //     rightmost_handle = _graph.get_handle(source);
    //     path_str = reverse_complement(path_str);
    //     path_spans_left_right = make_pair(path_spans_left_right.second, path_spans_left_right.first);
    // }
    // else
    // {
    //     leftmost_handle = _graph.get_handle(source);
    //     rightmost_handle = _graph.get_handle(sink);
    //     path_spans_left_right = path_spans_left_right;
    // }
    // cerr << "in move_path_to_snarl: " << endl;
    // cerr << "path_name: " << path_name << endl;
    // cerr << "path_str: " << path_str << endl;
    // cerr << "sink id and seq: " << sink << " " << _graph.get_sequence(_graph.get_handle(sink)) << endl;
    // cerr << "leftmost handle id and seq: " << _graph.get_id(leftmost_handle) << " " << _graph.get_sequence(leftmost_handle) << endl;


    // dealing with edge cases:
    // if source == sink, the snarl has become one single node and the path should map directly to that
    // one node.
    //todo!

    vector<pair<vector<handle_t>, int>> possible_path_starts = find_possible_path_starts(leftmost_handle, rightmost_handle, path_spans_left_right);
    // cerr << "size of possible_path_starts: " << possible_path_starts.size() << endl;
    vector<handle_t> new_path_location = extend_possible_paths(possible_path_starts, path_str, leftmost_handle, rightmost_handle, path_spans_left_right);
    // cerr << "size of new_path_location: " << new_path_location.size() << endl;


    // flip the order of the handles if the path moves right-to-left.
    if (!path_directed_left_to_right)
    {
        std::reverse(new_path_location.begin(), new_path_location.end());
        for (int i = 0; i != new_path_location.size(); i++)
        {
            // new_path_location[i] = _graph.app(new_path_location[i]);
            // cerr << "path_name: " << path_name << endl;
            // cerr << "handle_id: " << _graph.get_id(new_path_location[i]) << endl;
            // cerr << "handle seq: " << _graph.get_sequence(new_path_location[i]) << endl;
            // cerr << "is the handle reversed?: " << _graph.get_is_reverse(new_path_location[i]) << endl;
            // if (path_name != "CBS432.chrXIV")
            // {
                new_path_location[i] = _graph.flip(new_path_location[i]);
            // }
            // cerr << "is the handle reversed?: " << _graph.get_is_reverse(new_path_location[i]) << endl;
            // cerr << "is the handle at old_path in this pos reversed? " << _graph.get_is_reverse(old_path_location[i]) << endl;
        }
    }

    // cerr << "rewriting path " << _graph.get_path_name(_graph.get_path_handle_of_step(old_path.first)) << endl;
    // step_handle_t cur_old_step = old_path.first;
    // string old_path_series;
    // string old_path_str;
    // while (cur_old_step != old_path.second)
    // {
    //     cerr << "cur_old_step id: " << _graph.get_id(_graph.get_handle_of_step(cur_old_step)) << endl;
    //     old_path_series += " " + _graph.get_id(_graph.get_handle_of_step(cur_old_step));
    //     old_path_str += " " + _graph.get_sequence(_graph.get_handle_of_step(cur_old_step));
    //     // I wonder if I'm messing up the path sequence because get_sequence is always left->right, but the real path seq should be right->left on some handles.
    //     cur_old_step = _graph.get_next_step(cur_old_step);
    // }
    // string new_path_series;
    // string new_path_str;
    // for (handle_t handle : new_path_location)
    // {
    //     cerr << "cur_new_step id: " << _graph.get_id(handle) << endl;
    //     new_path_series += " " + _graph.get_id(handle);
    //     new_path_str += " " + _graph.get_sequence(handle);
    //     // I wonder if I'm messing up the path sequence because get_sequence is always left->right, but the real path seq should be right->left on some handles.
    // }

    // cerr << "request to rewrite segment with old_path: " << old_path_series << " " << old_path_str  << endl;
    // cerr << "and new path: " << new_path_series << " " << new_path_str  << endl;
    
    pair<step_handle_t, step_handle_t> new_path = _graph.rewrite_segment(old_path.first, old_path.second, new_path_location);

    // Test that the new path exists.
    if (new_path_location.size() == 0)
    {
        cerr << "************in UNIT_TEST for move_path_to_new_snarl************" << endl;
        cerr << "no new path location found." << endl;
    }
    // Test that the new path seq = old path seq.
    else
    {
        step_handle_t new_path_source;
        vector<step_handle_t> steps = _graph.steps_of_handle(new_path_location.front());
        for (auto step : steps) 
        {
            if (_graph.get_path_name(_graph.get_path_handle_of_step(step)) == path_name)
            {
                new_path_source = step;
            }
        }
        string new_path_str = _graph.get_sequence(_graph.get_handle_of_step(new_path_source));
        step_handle_t cur_step = new_path_source;
        while (_graph.get_id(_graph.get_handle_of_step(cur_step)) != _graph.get_id(new_path_location.back()))
        {
            //todo: make sure that this process correctly extracts path string.
            cur_step = _graph.get_next_step(cur_step);
            new_path_str += _graph.get_sequence(_graph.get_handle_of_step(cur_step)); 
        }
        string old_path_str;
        if (!path_directed_left_to_right)
        {
            old_path_str = reverse_complement(path_str);
        }
        else
        {
            old_path_str = path_str;    
        }
        if (old_path_str != new_path_str)
        {
            cerr << "************in UNIT_TEST for move_path_to_new_snarl************" << endl;
            cerr << "Once the path was moved into the new snarl, it didn't have the same sequence." << endl;
            cerr << "original seq: " << old_path_str << endl;
            cerr << "     new seq: " << new_path_str << endl;
        }
    }
    
    // cerr << "************END-UNIT_TEST for move_path_to_new_snarl.************"<< endl;

    return new_path;
}

// // Moves a path from its original location in the _graph to a new snarl,
// //      defined by a vector of interconnected handles.
// //      NOTE: the handles in new_snarl_handles may not preserve topological order after
// //      being passed to this method, if they were ordered before.
// // Arguments: _graph: the _graph containing the old_embedded_path and the handles in
// // new_snarl_topo_order
// //            old_embedded_path: a pair, where
// //                          pair.first is the first step_handle of interest in the
// //                          old_embedded_path, and pair.second is the step_handle *after*
// //                          the last step_handle of interest in the old_embedded_path (can
// //                          be the null step at the end of the path.)
// //            new_snarl_topo_order: all the handles in the new snarl, inside the _graph.
// // Return: None.
// void SnarlNormalizer::move_embedded_path_to_snarl(
//     const pair<step_handle_t, step_handle_t> &old_embedded_path,
//     vector<handle_t> &new_snarl_handles, id_t &new_source_id, id_t &new_sink_id,
//     const id_t &old_source_id, const id_t &old_sink_id, const bool backwards)
// {

//     // get the sequence associated with the path
//     string path_seq;
//     step_handle_t cur_step = old_embedded_path.first;

//     // if the old path is touching either/both the source/sink, we want to make sure that
//     // the newly moved path also touches those. Otherwise, any paths that extend beyond
//     // the source or sink may be cut into pieces when the portion of the path overlapping
//     // the snarl is moved to a region inside the snarl.
//     bool touching_source =
//         (_graph.get_id(_graph.get_handle_of_step(old_embedded_path.first)) ==
//             old_source_id);
//     bool touching_sink =
//         (_graph.get_id(_graph.get_handle_of_step(
//                 _graph.get_previous_step(old_embedded_path.second))) == old_sink_id);

//     cerr << "MOVE_PATH_TO_SNARL: touching_source and touching_sink: " << touching_source << " " << touching_sink << " new_source_id: " << new_source_id << " new_sink_id: " << new_sink_id << " old_source_id: " << old_source_id << " old_sink_id: " << old_sink_id << endl;
//     cerr << "old_embedded_path.second: " << _graph.get_id(_graph.get_handle_of_step(_graph.get_previous_step(old_embedded_path.second))) << endl;
//     // extract the path sequence of the embedded path:
//     while (cur_step != old_embedded_path.second)
//     {
//         path_seq += _graph.get_sequence(_graph.get_handle_of_step(cur_step));
//         cur_step = _graph.get_next_step(cur_step);
//     }
//     cerr << "path_seq to move: " << path_seq << endl;

//     vector<tuple<vector<handle_t>, int, int>> possible_paths;
//     for (handle_t handle : new_snarl_handles)
//     {
//         string handle_seq = _graph.get_sequence(handle);
//         if (backwards)
//         {
//             handle_seq = reverse_complement(handle_seq);
//         }

//         // starting index is where the path would begin in the handle,
//         // since it could begin in the middle of the handle.
//         vector<int> starting_indices =
//             check_handle_as_start_of_path_seq(handle_seq, path_seq);

//         cerr << "starting indices for path " << _graph.get_path_name(_graph.get_path_handle_of_step(old_embedded_path.first)) << "(path_seq: " << path_seq << ")"
//                 << " in handle " << _graph.get_id(handle) << " (handle_seq: " << handle_seq << ")" << endl;
//         for (int i : starting_indices)
//         {
//             cerr << i << endl;
//         }
//         cerr << "finished indices" << endl;

//         // if there is a starting index,
//         if (starting_indices.size() != 0)
//         {
//             for (int starting_index : starting_indices)
//             {
//                 if (((backwards && (starting_index >= path_seq.size())) || (!backwards && ((handle_seq.size() - starting_index) >= path_seq.size()))) &&
//                     source_and_sink_handles_map_properly(_graph, new_source_id,
//                                                             new_sink_id, touching_source,
//                                                             touching_sink, handle, handle))
//                 {
//                     // if the entire path fits inside the current handle, and if any
//                     // paths that touched source and sink in the old snarl would be
//                     // touching source and sink in the new snarl, then we've already
//                     // found the full mapping location of the path! Move the path, end
//                     // the method.
//                     vector<handle_t> new_path{handle};
//                     _graph.rewrite_segment(old_embedded_path.first,
//                                             old_embedded_path.second, new_path);
//                     // //todo: debug_statement
//                     // cerr << "found a full mapping at " << _graph.get_id(handle)
//                     //      << " w/ seq " << _graph.get_sequence(handle) << endl;
//                     return;
//                 }
//                 else
//                 {
//                     // this is a potential starting handle for the path. Add as a
//                     // possible_path.
//                     vector<handle_t> possible_path_handle_vec{handle};
//                     possible_paths.push_back(
//                         make_tuple(possible_path_handle_vec, starting_index,
//                                     handle_seq.size() - starting_index));
//                 }
//             }
//         }
//     }

//     // for every possible path, extend it to determine if it really is the path we're
//     // looking for:
//     while (!possible_paths.empty())
//     {
//         // take a path off of possible_paths, which will be copied for every iteration
//         // through _graph.follow_edges, below:
//         tuple<vector<handle_t>, int, int> possible_path_query = possible_paths.back();
//         possible_paths.pop_back();

//         // extend the path through all right-extending edges to see if any subsequent
//         // paths still satisfy the requirements for being a possible_path:
//         bool no_path = _graph.follow_edges(
//             get<0>(possible_path_query).back(), false, [&](const handle_t &next) {
//                 // make a copy to be extended for through each possible next handle in
//                 // follow edges.
//                 tuple<vector<handle_t>, int, int> possible_path = possible_path_query;

//                 // extract relevant information to make code more readable.
//                 string next_seq = _graph.get_sequence(next);
//                 id_t next_id = _graph.get_id(next);
//                 int &cur_index_in_path = get<2>(possible_path);
//                 if (cur_index_in_path <= path_seq.size() &&
//                     (find(new_snarl_handles.cbegin(), new_snarl_handles.cend(), next) !=
//                         new_snarl_handles.cend()))
//                 {
//                     // if the next handle would be the ending handle for the path,
//                     if (next_seq.size() >= (path_seq.size() - cur_index_in_path))
//                     {
//                         // cerr << "next handle would be the ending handle for the path"
//                         //      << endl;
//                         //     check to see if the sequence in the handle is suitable
//                         // for ending the path:
//                         int compare_length = path_seq.size() - cur_index_in_path;

//                         if ((next_seq.compare(0, compare_length, path_seq,
//                                                 cur_index_in_path, compare_length) == 0) &&
//                             source_and_sink_handles_map_properly(
//                                 _graph, new_source_id, new_sink_id, touching_source,
//                                 touching_sink, get<0>(possible_path).front(), next))
//                         {

//                             // we've found the new path! Move path to the new sequence,
//                             // and end the function.

//                             if (compare_length < next_seq.size())
//                             {
//                                 // If the path ends before the end of next_seq, then split
//                                 // the handle so that the path ends flush with the end of
//                                 // the first of the two split handles.

//                                 // divide the handle where the path ends;
//                                 pair<handle_t, handle_t> divided_next =
//                                     _graph.divide_handle(next, compare_length);
//                                 get<0>(possible_path).push_back(divided_next.first);

//                                 // Special case if next is the sink or the source, to
//                                 // preserve the reassignment of source and sink ids in
//                                 // integrate_snarl.
//                                 if (next_id == new_sink_id)
//                                 {
//                                     new_sink_id = _graph.get_id(divided_next.second);
//                                 }

//                                 // TODO: NOTE: finding the old "next" handle is expensive.
//                                 // TODO:   Use different container?
//                                 auto it = find(new_snarl_handles.begin(),
//                                                 new_snarl_handles.end(), next);

//                                 // replace the old invalidated handle with one of the new
//                                 // ones
//                                 *it = divided_next.first;
//                                 // stick the other new handle on the end of
//                                 // new_snarl_handles.
//                                 new_snarl_handles.push_back(divided_next.second);
//                             }
//                             else
//                             {
//                                 // otherwise, the end of the path already coincides with
//                                 // the end of the handle. In that case, just add it to the
//                                 // path.
//                                 get<0>(possible_path).push_back(next);
//                             }
//                             _graph.rewrite_segment(old_embedded_path.first,
//                                                     old_embedded_path.second,
//                                                     get<0>(possible_path));
//                             // //todo: debug_statement:
//                             // cerr << "got a full path: ";
//                             // for (handle_t handle : get<0>(possible_path)) {
//                             //     cerr << _graph.get_id(handle) << " ";
//                             // }
//                             // cerr << endl;

//                             // we've already found the path. No need to keep looking for
//                             // more paths.
//                             return false;
//                         }
//                     }
//                     // see if the next handle would be the continuation of the path, but
//                     // not the end,
//                     else
//                     {

//                         // check to see if the sequence in the handle is suitable for
//                         // extending the path:
//                         int compare_length = next_seq.size();
//                         // //todo: debug_statement
//                         // cerr << "compare returned false" << endl;
//                         // cerr << "compare in returned false: "
//                         //      << " next_seq len " << next_seq.size() << " compare_length
//                         //      "
//                         //      << compare_length << " path_seq len " << path_seq.size()
//                         //      << " cur_index_in_path " << cur_index_in_path << endl;
//                         // cerr << "if statement eval: cur_index_in_path <=
//                         // next_seq.size() "
//                         //      << (cur_index_in_path <= next_seq.size())
//                         //      << " next_seq.compare(0, compare_length, path_seq, "
//                         //         "cur_index_in_path, compare_length) == 0) "
//                         //      << (next_seq.compare(0, compare_length, path_seq,
//                         //                           cur_index_in_path, compare_length) ==
//                         //                           0)
//                         //      << endl;
//                         if (next_seq.compare(0, compare_length, path_seq,
//                                                 cur_index_in_path, compare_length) == 0)
//                         {
//                             // cerr << "compared in return false" << endl;
//                             // extend the path
//                             get<0>(possible_path).push_back(next);

//                             // update the current index in path_seq.
//                             get<2>(possible_path) += next_seq.size();

//                             // place back into possible_paths
//                             possible_paths.push_back(possible_path);
//                             // cerr << "extending the path!" << endl;
//                         }
//                     }
//                 }
//                 // continue to iterate through follow_edges.
//                 return true;
//             });

//         // //todo: debug_statement:
//         // if
//         // (graph.get_path_name(graph.get_path_handle_of_step(old_embedded_path.first))
//         // ==
//         //     "_alt_19f9bc9ad2826f58f113965edf36bb93740df46d_0") {
//         //     cerr << "mystery node 4214930: "
//         //          << _graph.get_sequence(graph.get_handle(4214930)) << endl;
//         // }

//         // if we've found a complete path in the above follow_edges, then we've
//         // already moved the path, and we're done.
//         if (!no_path)
//         {
//             return;
//         }
//     }
//     // //todo: figure out how to do some better error message instead of cerr.
//     // if we failed to find a path, show an error message.
//     cerr << "##########################\nWarning! Didn't find a corresponding path of "
//             "name "
//             << _graph.get_path_name(_graph.get_path_handle_of_step(old_embedded_path.first))
//             << " from the old snarl at " << old_source_id
//             << " in the newly aligned snarl. This snarl WILL be "
//             "normalized, resulting in a probably incorrectly-constructed snarl."
//             "\n##########################"
//             << endl
//             << endl;
//     // throw _graph.get_path_name(graph.get_path_handle_of_step(old_embedded_path.first));
//     // assert(true && "Warning! Didn't find a corresponding path of name " +
//     //         _graph.get_path_name(graph.get_path_handle_of_step(old_embedded_path.first))
//     //         + " from the old snarl in the newly aligned snarl.");
// }

} //VG
} //Algorithms
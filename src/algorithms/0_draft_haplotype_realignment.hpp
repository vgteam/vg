/*
Robin Rounthwaite
Find function call in ./subcommand/main.cpp
*/
#include "../gbwt_helper.hpp"
#include "../handle.hpp"
#include "../subgraph.hpp"
#include "../vg.hpp"
#include "count_walks.hpp"
#include <string>

/* TODO for improving haplotype_realignment.
Tomorrow:
* scale code upwards so that you can run code on every snarl in given graph.
* also add requirement that haps entering snarl = haps exiting snarl.
TODO: align haplotypes_not_at_source once we have a solution for alignments that insert
TODO:     the haplotype in a specified location
TODO:     (use more unique marker signals to identify where in other strings the
TODO:     middle-haplotype should align?)

TODO: consider splitting handles where embedded paths begin/end in the middle of a handle.
TODO:     (Note: would need to dynamically change other paths containing that handle. :-/)
TODO:     Or simply split the handles of interest and then realign the paths - expensive.
TODO:     Or insert *yet another* marker char to id where embedded paths begin/end, so its
TODO:     easily find where to split the handles afterwards. AND! it makes moving the
TODO:     paths less expensive.
TODO:     (fewer spots to check alignment in the snarl). If we have unique markers for
TODO:     each path, then
TODO:     it becomes O(N) time, instead of ~O(N*M*n) (N: number of bases in snarl; M:
TODO:     number of bases in path;
TODO:     n: number of potential starting places in the snarl (note: slightly less
TODO:     expensive since n is
TODO:     divided up among the M paths).)
TODO: this would also addres the possibility of an embedded path being moved to an
TODO:     alternative location
TODO:     when it overlaps a repetitive sequence. (previous thought, tho' above one is
TODO:     better): do I Need
TODO:      to account for this with a sense of "bases distant from source"?

TODO: make it so that gbwt file is customized by user rather than hardcoded.

TODO: make it so that you pass the gbwt file directory to a one-liner function
TODO:       (ran in normalize_main) that generates gbwt graph, extracts haps,
TODO:       aligns haps, and reintegrates haps. (eventually will do it for every
TODO:       snarl in the given graph).

*/
namespace vg {
void disambiguate_top_level_snarls(MutablePathDeletableHandleGraph &graph,
                                   const GBWTGraph &haploGraph, ifstream &snarl_stream);

bool disambiguate_snarl(MutablePathDeletableHandleGraph &graph,
                        const GBWTGraph &haploGraph, const id_t &source_id,
                        const id_t &sink_id);

pair<vector<vector<handle_t>>, vector<vector<handle_t>>>
extract_gbwt_haplotypes(const GBWTGraph &graph, const id_t &source_id,
                        const id_t &sink_id);

vector<vector<handle_t>>
find_haplotypes_not_at_source(const GBWTGraph &haploGraph,
                              unordered_set<handle_t> &touched_handles,
                              const id_t &sink_id);

vector<string> format_handle_haplotypes_to_strings(
    const GBWTGraph &haploGraph,
    const vector<vector<handle_t>> &haplotype_handle_vectors);

VG align_source_to_sink_haplotypes(const vector<string> &source_to_sink_haplotypes);

vector<pair<step_handle_t, step_handle_t>>
extract_embedded_paths_in_snarl(const PathHandleGraph &graph, const id_t &source_id,
                                const id_t &sink_id);

SubHandleGraph extract_subgraph(const HandleGraph &graph, const id_t &start_id,
                                const id_t &end_id);

void integrate_snarl(MutablePathDeletableHandleGraph &graph, const HandleGraph &new_snarl,
                     const vector<pair<step_handle_t, step_handle_t>> embedded_paths,
                     const id_t &source_id, const id_t &sink_id);

void move_path_to_snarl(MutablePathDeletableHandleGraph &graph,
                        const pair<step_handle_t, step_handle_t> &old_embedded_path,
                        vector<handle_t> &new_snarl_handles, id_t &source_id,
                        id_t &sink_id);

vector<int> check_handle_as_start_of_path_seq(const string &handle_seq,
                                              const string &path_seq);

// -------------------------------- DEBUG CODE BELOW: ------------------------------------

pair<vector<handle_t>, vector<handle_t>>
debug_get_sources_and_sinks(const HandleGraph &graph);

vector<string> debug_graph_to_strings(MutablePathDeletableHandleGraph &graph,
                                      id_t start_id, id_t end_id);

vector<string> debug_get_embedded_paths_from_source_to_sink(const PathHandleGraph &graph,
                                                            const handle_t &source_handle,
                                                            const handle_t &sink_handle);
} // namespace vg

/*
Deleted stuff:

void jordan_bug(MutablePathDeletableHandleGraph& graph){

    // example with one node:
    handle_t example = graph.get_handle(23448);
    handle_t replacement = graph.create_handle("GATTACA", 1);

    // move the source edges:
    //TODO: note the copy/paste. Ask if there's a better way to do this (I totally could
in Python!) graph.follow_edges(example, true,
                       [&](const handle_t &prev_handle) {
                           graph.create_edge(prev_handle, replacement);
                       });
    graph.follow_edges(example, false,
                       [&](const handle_t &next_handle) {
                           graph.create_edge(replacement, next_handle);
                       });

    // move the paths:
    graph.for_each_step_on_handle(example, [&](step_handle_t step)
    {
        graph.rewrite_segment(step, graph.get_next_step(step),
vector<handle_t>{replacement});
    });

    // example with two nodes:
    handle_t example_1 = graph.get_handle(23450);
    handle_t replacement_1 = graph.create_handle("GATTACA", 2);
    handle_t replacement_2 = graph.create_handle("GATTACA", 3);
    graph.create_edge(replacement_1, replacement_2);

    // move the source edges:
    //TODO: note the copy/paste. Ask if there's a better way to do this (I totally could
in Python!) graph.follow_edges(example_1, true,
                       [&](const handle_t &prev_handle) {
                           graph.create_edge(prev_handle, replacement_1);
                       });
    graph.follow_edges(example_1, false,
                       [&](const handle_t &next_handle) {
                           graph.create_edge(replacement_2, next_handle);
                       });

    // move the paths:
    graph.for_each_step_on_handle(example_1, [&](step_handle_t step)
    {
        graph.rewrite_segment(step, step, vector<handle_t>{replacement_1, replacement_2});
    });
}
 */

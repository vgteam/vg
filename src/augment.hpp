#ifndef VG_AUGMENT_HPP_INCLUDED
#define VG_AUGMENT_HPP_INCLUDED

#include <iostream>
#include <sstream>
#include <string>
#include <map>
#include <limits>
#include <functional>

#include "handle.hpp"

namespace vg {

class Packer;

using namespace std;

/// %Edit the graph to include all the sequence and edges added by the given
/// paths. Can handle paths that visit nodes in any orientation. Note that
/// this method sorts the graph and rebuilds the path index, so it should
/// not be called in a loop.
///
/// if gam_path is "-", then stdin used
/// if gam_out_path is "-", then stdout used
/// If gam_out_path is not empty, the paths will be modified to reflect their
/// embedding in the modified graph and written to the path.
/// aln_format used to toggle between GAM and GAF
/// If out_translation is not null, a list of translations, one per node existing
/// after the edit, describing
/// how each new or conserved node is embedded in the old graph. 
/// if embed_paths is true, then the augmented alignemnents will be saved as embededed paths in the graph
/// in order to add it back to the graph.
/// If break_at_ends is true, nodes will be broken at
/// the ends of paths that start/end woth perfect matches, so the paths can
/// be added to the vg graph's paths object.
/// If soft_clip is true, soft clips will be removed from the input paths
/// before processing, and the dangling ends won't end up in the graph
/// If filter_out_of_graph_alignments is true, some extra time will be taken to check if
/// all nodes in the alignment are in the graph.  If they aren't, then it will be ignored
/// If an edit sequence's avg base quality is less than min_baseq it will be ignored (considered a match)
/// If an alignment's mapping quality is less than min_mapq it is ignored
/// A packer is required for all non-mapq filters
/// If a breakpoint has less than min_bp_coverage it is not included in the graph
/// Edits with more than max_frac_n N content will be ignored
void augment(MutablePathMutableHandleGraph* graph,
             const string& gam_path,
             const string& aln_format = "GAM",
             vector<Translation>* out_translation = nullptr,
             const string& gam_out_path = "",
             bool embed_paths = false,
             bool break_at_ends = false,
             bool remove_soft_clips = false,
             bool filter_out_of_graph_alignments = false,
             double min_baseq = 0,
             double min_mapq = 0,
             Packer* packer = nullptr,
             size_t min_bp_coverage = 0,
             double max_frac_n = 1.,
             bool edges_only = false);

/// Like above, but operates on a vector of Alignments, instead of a file
/// (Note: It is best to use file interface to stream large numbers of alignments to save memory)
void augment(MutablePathMutableHandleGraph* graph,
             vector<Path>& path_vector,
             const string& aln_format = "GAM",
             vector<Translation>* out_translation = nullptr,
             const string& gam_out_path = "",
             bool embed_paths = false,
             bool break_at_ends = false,
             bool remove_soft_clips = false,
             bool filter_out_of_graph_alignments = false,
             double min_baseq = 0,
             double min_mapq = 0,
             Packer* packer = nullptr,
             size_t min_bp_coverage = 0,
             double max_frac_n = 1.,
             bool edges_only = false);

/// Generic version used to implement the above three methods.  
void augment_impl(MutablePathMutableHandleGraph* graph,
                  function<void(function<void(Alignment&)>, bool, bool)> iterate_gam,
                  const string& aln_format,
                  vector<Translation>* out_translation,
                  const string& gam_out_path,
                  bool embed_paths,
                  bool break_at_ends,
                  bool remove_soft_clips,
                  bool filter_out_of_graph_alignments,
                  double min_baseq,
                  double min_mapq,
                  Packer* packer,
                  size_t min_bp_coverage,
                  double max_frac_n,
                  bool edges_only);

/// Add a path to the graph.  This is like VG::extend, and expects
/// a path with no edits, and for all the nodes and edges in the path
/// to exist exactly in the graph
path_handle_t add_path_to_graph(MutablePathHandleGraph* graph, const Path& path);

/// Compute the average base quality of an edit.
/// If the edit has no sequence or there are no base_quals given,
/// then double_max is returned. 
double get_avg_baseq(const Edit& edit, const string& base_quals, size_t position_in_read);

/// Find all the points at which a Path enters or leaves nodes in the graph. Adds
/// them to the given map by node ID of sets of bases in the node that will need
/// to become the starts of new nodes.
///
/// If break_ends is true, emits breakpoints at the ends of the path, even
/// if it starts/ends with perfect matches.

/// Find all the points at which a Path enters or leaves nodes in the graph. Adds
/// them to the given map by node ID of sets of bases in the node that will need
/// to become the starts of new nodes.
///
/// If break_ends is true, emits breakpoints at the ends of the path, even
/// if it starts/ends with perfect matches.
void find_breakpoints(const Path& path, unordered_map<id_t, set<pos_t>>& breakpoints, bool break_ends = true,
                      const string& base_quals = "", double min_baseq = 0, double max_frac_n = 1.);

/// Flips the breakpoints onto the forward strand.
unordered_map<id_t, set<pos_t>> forwardize_breakpoints(const HandleGraph* graph,
                                                       const unordered_map<id_t, set<pos_t>>& breakpoints);


/// Like "find_breakpoints", but store in packed structure (better for large gams and enables coverage filter)
void find_packed_breakpoints(const Path& path, Packer& packed_breakpoints, bool break_ends = true,
                             const string& base_quals = "", double min_baseq = 0, double max_frac_n = 1.);

/// Filters the breakpoints by coverage, and converts them back from the Packer to the STL map
/// expected by following methods
unordered_map<id_t, set<pos_t>> filter_breakpoints_by_coverage(const Packer& packed_breakpoints, size_t min_bp_coverage);

/// Take a map from node ID to a set of offsets at which new nodes should
/// start (which may include 0 and 1-past-the-end, which should be ignored),
/// break the specified nodes at those positions. Returns a map from old
/// node start position to new node pointer in the graph. Note that the
/// caller will have to crear and rebuild path rank data.
///
/// Returns a map from old node start position to new node. This map
/// contains some entries pointing to null, for positions past the ends of
/// original nodes. It also maps from positions on either strand of the old
/// node to the same new node pointer; the new node's forward strand is
/// always the same as the old node's forward strand.
map<pos_t, id_t> ensure_breakpoints(MutableHandleGraph* graph,
                                    const unordered_map<id_t, set<pos_t>>& breakpoints);

/// Remove edits in our graph that don't correspond to breakpoints (ie were effectively filtered
/// out due to insufficient coverage.  This way, subsequent logic in add_nodes_and_edges
/// can be run correctly.  Returns true if at least one edit survived the filter.
bool simplify_filtered_edits(HandleGraph* graph, Alignment& aln, Path& path, const map<pos_t, id_t>& node_translation,
                             const unordered_map<id_t, size_t>& orig_node_sizes,
                             double min_baseq = 0, double max_frac_n = 1.);

/// Given a path on nodes that may or may not exist, and a map from start
/// position in the old graph to a node in the current graph, add all the
/// new sequence and edges required by the path. The given path must not
/// contain adjacent perfect match edits in the same mapping, or any
/// deletions on the start or end of mappings (the removal of which can be
/// accomplished with the Path::simplify() function).
///
/// Outputs (and caches for subsequent calls) novel node runs in added_seqs,
/// and Paths describing where novel nodes translate back to in the original
/// graph in added_nodes. Also needs a map of the original sizes of nodes
/// deleted from the original graph, for reverse complementing. If dangling
/// is nonempty, left edges of nodes created for initial inserts will
/// connect to the specified sides. At the end, dangling is populated with
/// the side corresponding to the last edit in the path.
///
/// Returns a fully embedded version of the path, after all node insertions,
/// divisions, and translations.
Path add_nodes_and_edges(MutableHandleGraph* graph,
                         const Path& path,
                         const map<pos_t, id_t>& node_translation,
                         unordered_map<pair<pos_t, string>, vector<id_t>>& added_seqs,
                         unordered_map<id_t, Path>& added_nodes,
                         const unordered_map<id_t, size_t>& orig_node_sizes,
                         set<NodeSide>& dangling,
                         size_t max_node_size = 1024);
    
/// This version doesn't require a set of dangling sides to populate                         
Path add_nodes_and_edges(MutableHandleGraph* graph,
                         const Path& path,
                         const map<pos_t, id_t>& node_translation,
                         unordered_map<pair<pos_t, string>, vector<id_t>>& added_seqs,
                         unordered_map<id_t, Path>& added_nodes,
                         const unordered_map<id_t, size_t>& orig_node_sizes,
                         size_t max_node_size = 1024);

/// Produce a graph Translation object from information about the editing process.
vector<Translation> make_translation(const HandleGraph* graph,
                                     const map<pos_t, id_t>& node_translation,
                                     const unordered_map<id_t, Path>& added_nodes,
                                     const unordered_map<id_t, size_t>& orig_node_sizes);
  
/// Add edges between consecutive mappings that aren't already in the graph
/// note: offsets are completely ignored (a simplifying assumption designed
/// to help with SV genotpying with pack/call as edge packing works similarly)
///
/// No existing nodes or edges are modified, and no nodes are added, just edges
/// So the output graph will be id-space compatible, and any GAM/GAF will continue
/// to be valid for it.  
void add_edges_only(MutableHandleGraph* graph,
                    function<void(function<void(Alignment&)>,bool, bool)> iterate_gam,
                    double min_mapq,
                    size_t min_bp_coverage);

}

#endif

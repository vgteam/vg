#ifndef VG_PATH_HPP_INCLUDED
#define VG_PATH_HPP_INCLUDED

#include <iostream>
#include <algorithm>
#include <functional>
#include <set>
#include <list>
#include <sstream>
#include "json2pb.h"
#include "vg.pb.h"
#include "edit.hpp"
#include "hash_map.hpp"
#include "utility.hpp"
#include "types.hpp"
#include "position.hpp"
#include "nodetraversal.hpp"

//#define debug

namespace vg {

using namespace std;

class Paths {
public:

    Paths(void);

    // copy
    Paths(const Paths& other) {
        if (this != &other) {
            _paths = other._paths;
            rebuild_node_mapping();
        }
    }
    // move constructor
    Paths(Paths&& other) noexcept {
        _paths = other._paths;
        other.clear();
        rebuild_node_mapping();
    }

    // copy assignment operator
    Paths& operator=(const Paths& other) {
        Paths tmp(other);
        *this = std::move(tmp);
        return *this;
    }

    // move assignment operator
    Paths& operator=(Paths&& other) noexcept {
        std::swap(_paths, other._paths);
        rebuild_node_mapping();
        return *this;
    }

    // This maps from path name to the list of Mappings for that path.
    map<string, list<Mapping> > _paths;
    // This maps from Mapping* pointer to its iterator in its list of Mappings
    // for its path. The list in question is stored above in _paths. Recall that
    // std::list iterators are bidirectional.
    map<Mapping*, list<Mapping>::iterator > mapping_itr;
    // This maps from Mapping* pointer to the name of the path it belongs to
    // (which can then be used to get the list its iterator belongs to).
    map<Mapping*, string> mapping_path;
    void sort_by_mapping_rank(void);
    /// Reassign ranks and rebuild indexes, treating the mapping lists in _paths as the truth.
    void rebuild_mapping_aux(void);
    // We need this in order to make sure we aren't adding duplicate mappings
    // with the same rank in the same path. Maps from path name and rank to
    // Mapping pointer.
    map<string, map<size_t, Mapping*>> mappings_by_rank;
    // This maps from node ID, then path name, then rank and orientation, to
    // Mapping pointers for the mappings on that path to that node.
    map<id_t, map<string, set<Mapping*>>> node_mapping;
    // record which head nodes we have
    // we'll use this when determining path edge crossings--- all paths implicitly cross these nodes
    set<id_t> head_tail_nodes;
    bool is_head_or_tail_node(id_t);
    vector<string> all_path_names(void);
    // records which paths are circular
    set<string> circular;
    void make_circular(const string& name);
    void make_linear(const string& name);
    
    void rebuild_node_mapping(void);
    
    // Find the given mapping in its path, so mappings can be inserted before
    // it.
    list<Mapping>::iterator find_mapping(Mapping* m);
    // remove the given Mapping from its path. Returns an iterator to the
    // mapping that came after it, or the end of the list if no mapping came
    // after it.
    list<Mapping>::iterator remove_mapping(Mapping* m);
    // insert the given mapping into the given path, before the mapping pointed
    // to by the given iterator. Returns an iterator to the newly-inserted
    // mapping.
    list<Mapping>::iterator insert_mapping(list<Mapping>::iterator w,
                                           const string& path_name, const Mapping& m);
    pair<Mapping*, Mapping*> divide_mapping(Mapping* m, const Position& pos);
    pair<Mapping*, Mapping*> divide_mapping(Mapping* m, size_t offset);
    // replace the mapping with two others in the order provided
    pair<Mapping*, Mapping*> replace_mapping(Mapping* m, pair<Mapping, Mapping> n);
    // Note that this clears and rebuilds all the indexes
    void remove_paths(const set<string>& names);
    // This one actually unthreads the path from the indexes. It's O(path
    // length), but you can do several calls without thrashing the index
    // clear and rebuild functions.
    void remove_path(const string& name);
    void keep_paths(const set<string>& name);
    void remove_node(id_t id);
    bool has_path(const string& name);
    void to_json(ostream& out);
    list<Mapping>& get_path(const string& name);
    list<Mapping>& get_create_path(const string& name);
    list<Mapping>& create_path(const string& name);
    // Does the given path have a mapping meeting the given criteria?
    // Is there a mapping in the given path with the given assigned rank? Note
    // that the rank passed may not be 0.
    bool has_mapping(const string& name, size_t rank);
    // We used to be able to search for a Mapping by value, but that's not
    // efficient if the Mappings don't have ranks, and it never checked the
    // edits for equality anyway.
    bool has_node_mapping(id_t id);
    bool has_node_mapping(Node* n);
    map<string, set<Mapping*> >& get_node_mapping(Node* n);
    map<string, set<Mapping*> >& get_node_mapping(id_t id);
    map<string, map<int, Mapping*> > get_node_mappings_by_rank(id_t id);
    // Copy all the mappings for the node with the given ID, and return them in
    // a map by path name and then by rank.
    map<string, map<int, Mapping> > get_node_mapping_copies_by_rank(id_t id);
    // Go left along the path that this Mapping* belongs to, and return the
    // Mapping* there, or null if this Mapping* is the first in its path.
    Mapping* traverse_left(Mapping* mapping);
    // Go right along the path that this Mapping* belongs to, and return the
    // Mapping* there, or null if this Mapping* is the last in its path.
    Mapping* traverse_right(Mapping* mapping);
    // TODO: should this be a reference?
    const string mapping_path_name(Mapping* m);
    // the patsh of the node
    set<string> of_node(id_t id);
    // get the paths on this node and the number of mappings from each one
    map<string, int> node_path_traversal_counts(id_t id, bool rev = false);
    // the return vector contains the name of each path that crosses the node
    // repeated as many times as it crosses
    vector<string> node_path_traversals(id_t id, bool rev = false);
    bool are_consecutive_nodes_in_path(id_t id1, id_t id2, const string& path_name);
    // paths that cross the edge, 
    vector<string> over_edge(id_t id1, bool rev1,
                             id_t id2, bool rev2,
                             vector<string> following);
    vector<string> over_directed_edge(id_t id1, bool rev1,
                                      id_t id2, bool rev2,
                                      vector<string> following);
    size_t size(void) const;
    bool empty(void) const;
    // clear the internal data structures tracking mappings and storing the paths
    void clear(void);
    void clear_mapping_ranks(void);
    void compact_ranks(void);
    //void add_node_mapping(Node* n);
    void load(istream& in);
    void write(ostream& out);
    /// Add all paths into the given Protobuf graph. Creates a new path for every path.
    void to_graph(Graph& g);
    // get a path
    Path path(const string& name);
    // add mappings. Mappings are assumed to either be pre-sorted or to have
    // ranks assigned already. If you append a mapping with no rank after
    // mappings with ranks, or visa versa, the relative ordering is undefined.
    // Also, if both pre-ranked and un-ranked mappings are appended, or if you
    // called clear_mapping_ranks before appending, mappings_by_rank may be
    // incorrect until rebuild_mapping_aux() is called.
    void append_mapping(const string& name, const Mapping& m);
    void append_mapping(const string& name, id_t id, size_t rank = 0, bool is_reverse = false);
    void prepend_mapping(const string& name, const Mapping& m);
    void prepend_mapping(const string& name, id_t id, size_t rank = 0, bool is_reverse = false);
    size_t get_next_rank(const string& name);
    void append(Paths& p);
    void append(Graph& g);
    void extend(Paths& p);
    void extend(const Path& p);
    void for_each(const function<void(const Path&)>& lambda);
    // Loop over the names of paths without actually extracting the Path objects.
    void for_each_name(const function<void(const string&)>& lambda);
    void for_each_stream(istream& in, const function<void(Path&)>& lambda);
    void increment_node_ids(id_t inc);
    // Replace the node IDs used as keys with those used as values.
    // This is only efficient to do in a batch.
    void swap_node_ids(hash_map<id_t, id_t>& id_mapping);
    // sets the mapping to the new id
    // erases current (old index information)
    void reassign_node(id_t new_id, Mapping* m);
    void for_each_mapping(const function<void(Mapping*)>& lambda);
};

string  path_to_string(Path p);
Path& increment_node_mapping_ids(Path& p, id_t inc);
Path& append_path(Path& a, const Path& b);
const Paths paths_from_graph(Graph& g);
void parse_region(const string& target, string& name, id_t& start, id_t& end);
int path_to_length(const Path& path);
int path_from_length(const Path& path);
int mapping_to_length(const Mapping& m);
int mapping_from_length(const Mapping& m);
Position first_path_position(const Path& path);
Position last_path_position(const Path& path);
int to_length(const Mapping& m);
int from_length(const Mapping& m);
bool mapping_ends_in_deletion(const Mapping& m);
bool mapping_starts_in_deletion(const Mapping& m);
bool mapping_is_total_deletion(const Mapping& m);
bool mapping_is_simple_match(const Mapping& m);
bool path_is_simple_match(const Path& p);
// convert the mapping to the particular node into the sequence implied by the mapping
const string mapping_sequence(const Mapping& m, const string& node_seq);
const string mapping_sequence(const Mapping& m, const Node& n);
// Reverse-complement a Mapping and all the Edits in it. A function to get node
// lengths is needed, because the mapping will need to count its position from
// the other end of the node.
Mapping reverse_complement_mapping(const Mapping& m,
                                   const function<int64_t(id_t)>& node_length);
// Reverse-complement a Path and all the Mappings in it. A function to get node
// lengths is needed, because the mappings will need to count their positions
// from the other ends of their nodes.
Path reverse_complement_path(const Path& path,
                             const function<int64_t(id_t)>& node_length);
/// Simplify the path for addition as new material in the graph. Remove any
/// mappings that are merely single deletions, merge adjacent edits of the same
/// type, strip leading and trailing deletion edits on mappings, and make sure no
/// mappings have missing positions.
Path simplify(const Path& p);
/// Merge adjacent edits of the same type, strip leading and trailing deletion
/// edits (while updating positions if necessary), and makes sure position is
/// actually set.
Mapping simplify(const Mapping& m);

/// Merge adjacent edits of the same type
Path merge_adjacent_edits(const Path& m);
/// Merge adjacent edits of the same type
Mapping merge_adjacent_edits(const Mapping& m);
// trim path so it starts and begins with a match (or is empty)
Path trim_hanging_ends(const Path& p);
// make a new mapping that concatenates the mappings
Mapping concat_mappings(const Mapping& m, const Mapping& n);
// make a new path that concatenates the two given paths
Path concat_paths(const Path& path1, const Path& path2);
// extend the first path by the second, avoiding copy operations
Path& extend_path(Path& path1, const Path& path2);
// divide mapping at reference-relative position
pair<Mapping, Mapping> cut_mapping(const Mapping& m, const Position& pos);
// divide mapping at target-relative offset (as measured in to_length)
pair<Mapping, Mapping> cut_mapping(const Mapping& m, size_t offset);
// divide path at reference-relative position
pair<Path, Path> cut_path(const Path& path, const Position& pos);
// divide the path at a path-relative offset as measured in to_length from start
pair<Path, Path> cut_path(const Path& path, size_t offset);
bool maps_to_node(const Path& p, id_t id);
// the position that starts just after the path ends
Position path_start(const Path& path);
Position path_end(const Path& path);
bool adjacent_mappings(const Mapping& m1, const Mapping& m2);
// Return true if a mapping is a perfect match (i.e. contains no non-match edits)
bool mapping_is_match(const Mapping& m);
double divergence(const Mapping& m);
// Return the identity for the path: perfect matches over total length.
// For zero-length paths, returns 0.
double identity(const Path& path);
// compare the agreement between two alignments
double overlap(const Path& p1, const Path& p2);
// helps estimate overapls quickly
void decompose(const Path& path, map<pos_t, int>& ref_positions, map<pos_t, Edit>& edits);

// switches the node ids in the path to the ones indicated by the translator
void translate_node_ids(Path& path, const unordered_map<id_t, id_t>& translator);
// switches the node ids and orientations in the path to the ones indicated by the translator
void translate_oriented_node_ids(Path& path, const unordered_map<id_t, pair<id_t, bool>>& translator);
    
// the first position on the path
pos_t initial_position(const Path& path);
// the last position on the path
pos_t final_position(const Path& path);
    
// Turn a list of node traversals into a path
Path path_from_node_traversals(const list<NodeTraversal>& traversals);

}

#endif

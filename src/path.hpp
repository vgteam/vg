#ifndef VG_PATH_HPP_INCLUDED
#define VG_PATH_HPP_INCLUDED

#include <iostream>
#include <algorithm>
#include <functional>
#include <set>
#include <list>
#include <sstream>
#include "vg/io/json2pb.h"
#include <vg/vg.pb.h>
#include "vg/io/edit.hpp"
#include "hash_map.hpp"
#include "utility.hpp"
#include "types.hpp"
#include "position.hpp"
#include "nodetraversal.hpp"

//#define debug

namespace vg {

using namespace std;

class mapping_t {
public:
    int64_t traversal; // negative implies reverse
    int32_t length;
    int32_t rank;
    mapping_t(void);
    mapping_t(const Mapping& m);
    Mapping to_mapping(void) const;
    id_t node_id(void) const;
    void set_node_id(id_t id);
    bool is_reverse(void) const;
    void set_is_reverse(bool is_rev);
};

/// Allow a mapping_t to be printed, for debugging purposes
ostream& operator<<(ostream& out, mapping_t mapping);

class Paths {
public:

    // This predicate matches the names of alt paths.
    // We used to use a regex but that's a very slow way to check a prefix.
    const static function<bool(const string&)> is_alt;

    // Use the path metadata api to strip a subrange from a path name.
    // If desired, the subrange can be stored in start and end
    static string strip_subrange(const string& path_name, subrange_t* subrange = nullptr);

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
        _paths = std::move(other._paths);
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
    map<string, list<mapping_t> > _paths;

private:
    // These need to be private because path names are lazily assigned IDs by the accessors
    // They also need to be mutable because we want our accessors to treat our object as logical const
    mutable int64_t max_path_id;
    mutable map<string, int64_t> name_to_id;
    mutable map<int64_t, string> id_to_name;
public:
    /// Get the lazily assigned numeric ID for a path, by name.
    int64_t get_path_id(const string& name) const;
    /// Get the name of a path, by numeric ID.
    const string& get_path_name(int64_t id) const;
    
    // This maps from mapping_t* pointer to its iterator in its list of Mappings
    // for its path and the id of the path.
    // The list in question is stored above in _paths.
    // Recall that std::list iterators are bidirectional.
    hash_map<const mapping_t*, pair<list<mapping_t>::iterator, int64_t> > mapping_itr;
    // This maps from mapping_t* pointer to the name of the path it belongs to
    // (which can then be used to get the list its iterator belongs to).
    void sort_by_mapping_rank(void);
    /// Reassign ranks and rebuild indexes, treating the mapping lists in _paths as the truth.
    void rebuild_mapping_aux(void);
    // We need this in order to make sure we aren't adding duplicate mappings
    // with the same rank in the same path. Maps from path name and rank to
    // Mapping pointer.
    map<string, hash_map<int32_t, mapping_t*>> mappings_by_rank;
    // This maps from node ID, then path name, then rank and orientation, to
    // Mapping pointers for the mappings on that path to that node.
    hash_map<id_t, map<int64_t, set<mapping_t*>>> node_mapping;
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
    list<mapping_t>::iterator find_mapping(mapping_t* m);
    // remove the given Mapping from its path. Returns an iterator to the
    // mapping that came after it, or the end of the list if no mapping came
    // after it.
    list<mapping_t>::iterator remove_mapping(mapping_t* m);
    // insert the given mapping into the given path, before the mapping pointed
    // to by the given iterator. Returns an iterator to the newly-inserted
    // mapping.
    list<mapping_t>::iterator insert_mapping(list<mapping_t>::iterator w,
                                             const string& path_name, const mapping_t& m);
    pair<mapping_t*, mapping_t*> divide_mapping(mapping_t* m, const Position& pos);
    pair<mapping_t*, mapping_t*> divide_mapping(mapping_t* m, size_t offset);
    // replace the mapping with two others in the order provided
    pair<mapping_t*, mapping_t*> replace_mapping(mapping_t* m, pair<mapping_t, mapping_t> n);
    // Note that this clears and rebuilds all the indexes
    void remove_paths(const set<string>& names);
    // This one actually unthreads the path from the indexes. It's O(path
    // length), but you can do several calls without thrashing the index
    // clear and rebuild functions.
    void remove_path(const string& name);
    void keep_paths(const set<string>& name);
    void remove_node(id_t id);
    bool has_path(const string& name) const;
    void to_json(ostream& out);
    list<mapping_t>& get_path(const string& name);
    list<mapping_t>& get_create_path(const string& name);
    list<mapping_t>& create_path(const string& name);
    // Does the given path have a mapping meeting the given criteria?
    // Is there a mapping in the given path with the given assigned rank? Note
    // that the rank passed may not be 0.
    bool has_mapping(const string& name, int32_t rank);
    // We used to be able to search for a Mapping by value, but that's not
    // efficient if the Mappings don't have ranks, and it never checked the
    // edits for equality anyway.
    bool has_node_mapping(id_t id) const;
    bool has_node_mapping(Node* n) const;
    map<int64_t, set<mapping_t*> >& get_node_mapping(Node* n);
    map<int64_t, set<mapping_t*> >& get_node_mapping(id_t id);
    const map<int64_t, set<mapping_t*> >& get_node_mapping(id_t id) const;
    map<string, set<mapping_t*> > get_node_mapping_by_path_name(Node* n);
    map<string, set<mapping_t*> > get_node_mapping_by_path_name(id_t id);
    map<string, map<int, mapping_t*> > get_node_mappings_by_rank(id_t id);
    // Copy all the mappings for the node with the given ID, and return them in
    // a map by path name and then by rank.
    map<string, map<int, mapping_t> > get_node_mapping_copies_by_rank(id_t id);
    // Go left along the path that this mapping_t* belongs to, and return the
    // mapping_t* there, or null if this Mapping* is the first in its path.
    mapping_t* traverse_left(mapping_t* mapping);
    // Go right along the path that this mapping_t* belongs to, and return the
    // mapping_t* there, or null if this mapping_t* is the last in its path.
    mapping_t* traverse_right(mapping_t* mapping);
    // TODO: should this be a reference?
    const string mapping_path_name(mapping_t* m);
    int64_t mapping_path_id(mapping_t* m);
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
    /// Write paths to the given stream, ending with an EOF marker.
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
    void append_mapping(const string& name, const mapping_t& m, bool warn_on_duplicates = false);
    void append_mapping(const string& name, id_t id, bool is_reverse, size_t length, size_t rank = 0, bool warn_on_duplicates = false);
    // TODO: Adapt this to use mapping_t instead.
    void prepend_mapping(const string& name, const Mapping& m, bool warn_on_duplicates = false);
    void prepend_mapping(const string& name, id_t id, bool is_reverse, size_t length, size_t rank = 0, bool warn_on_duplicates = false);
    size_t get_next_rank(const string& name);
    void append(const Paths& paths, bool warn_on_duplicates = false, bool rebuild_indexes = true);
    void append(const Graph& g, bool warn_on_duplicates = false, bool rebuild_indexes = true);
    void extend(const Paths& paths, bool warn_on_duplicates = false, bool rebuild_indexes = true);
    void extend(const Path& p, bool warn_on_duplicates = false, bool rebuild_indexes = true);
    void extend(const vector<Path> & paths, bool warn_on_duplicates = false, bool rebuild_indexes = true);
    void for_each(const function<void(const Path&)>& lambda);
    // Loop over the names of paths without actually extracting the Path objects.
    void for_each_name(const function<void(const string&)>& lambda) const;
    // Like for_each_name but allows stopping early.
    // TODO: Use the libhandlegraph unified iteratee pattern here.
    bool for_each_name_stoppable(const function<bool(const string&)>& lambda) const;
    void for_each_stream(istream& in, const function<void(Path&)>& lambda);
    void increment_node_ids(id_t inc);
    // Replace the node IDs according to a mapping from old ID to new ID.
    // This is only efficient to do in a batch.
    void swap_node_ids(const std::function<nid_t(const nid_t&)>& get_new_id);
    void swap_node_ids(hash_map<id_t, id_t>& id_mapping);
    // sets the mapping to the new id
    // erases current (old index information)
    void reassign_node(id_t new_id, mapping_t* m);
    void for_each_mapping(const function<void(mapping_t&)>& lambda);
};

string  path_to_string(Path p);
Path& increment_node_mapping_ids(Path& p, id_t inc);
/// Append the second path onto the end of the first, without combining
/// mappings or simplifying. Modifies and returns a reference to the first
/// path.
Path& append_path(Path& a, const Path& b);
const Paths paths_from_graph(Graph& g);
int path_to_length(const Path& path);
int path_from_length(const Path& path);
int mapping_to_length(const Mapping& m);
int mapping_from_length(const Mapping& m);
Position first_path_position(const Path& path);
Position last_path_position(const Path& path);
int to_length(const Mapping& m);
int from_length(const Mapping& m);
bool mappings_equivalent(const Mapping& m1, const Mapping& m2);
bool mapping_ends_in_deletion(const Mapping& m);
bool mapping_starts_in_deletion(const Mapping& m);
bool mapping_is_total_insertion(const Mapping& m);
bool mapping_is_total_deletion(const Mapping& m);
bool mapping_is_simple_match(const Mapping& m);
bool path_is_simple_match(const Path& p);
// convert the mapping to the particular node into the sequence implied by the mapping
const string mapping_sequence(const Mapping& m, const string& node_seq);
const string mapping_sequence(const Mapping& m, const Node& n);
// convert the path to a sequence
string path_sequence(const HandleGraph& graph, const Path& path);
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
// Reverse-complement a Mapping and all the Edits in it in place.
void reverse_complement_mapping_in_place(Mapping* m,
                                         const function<int64_t(id_t)>& node_length);
// Reverse-complement a Path and all the Mappings in it in place.
void reverse_complement_path_in_place(Path* path,
                                      const function<int64_t(id_t)>& node_length);
/// Simplify the path for addition as new material in the graph. Remove any
/// mappings that are merely single deletions, merge adjacent edits of the same
/// type, strip leading and trailing deletion edits on mappings, and make sure no
/// mappings have missing positions.
Path simplify(const Path& p, bool trim_internal_deletions = true);
/// Merge adjacent edits of the same type, strip leading and trailing deletion
/// edits (while updating positions if necessary), and makes sure position is
/// actually set.
Mapping simplify(const Mapping& m, bool trim_internal_deletions = true);

/// Return true if two edits could be combined into one (assuming adjacency).
bool edits_are_compatible(const Edit& e, const Edit& f);
/// Glom the second edit into the first, assuming adjacency.
void merge_edits_in_place(Edit& e, const Edit& f);

/// Merge adjacent edits of the same type
Path merge_adjacent_edits(const Path& m);
/// Merge adjacent edits of the same type
Mapping merge_adjacent_edits(const Mapping& m);
// trim path so it starts and begins with a match (or is empty)
Path trim_hanging_ends(const Path& p);
// make a new mapping that concatenates the mappings
Mapping concat_mappings(const Mapping& m, const Mapping& n, bool trim_internal_deletions = true);
// make a new path that concatenates the two given paths
Path concat_paths(const Path& path1, const Path& path2);
// extend the first path by the second, avoiding copy operations
Path& extend_path(Path& path1, const Path& path2);
// divide mapping at reference-relative position
pair<Mapping, Mapping> cut_mapping(const Mapping& m, const Position& pos);
pair<mapping_t, mapping_t> cut_mapping(const mapping_t& m, const Position& pos);
// divide mapping at reference-relative offset (as measure in from_length)
pair<Mapping, Mapping> cut_mapping_offset(const Mapping& m, size_t offset);
pair<mapping_t, mapping_t> cut_mapping_offset(const mapping_t& m, size_t offset);
// divide mapping at target-relative offset (as measured in to_length)
pair<Mapping, Mapping> cut_mapping(const Mapping& m, size_t offset);
pair<mapping_t, mapping_t> cut_mapping(const mapping_t& m, size_t offset);
// divide path at reference-relative position
pair<Path, Path> cut_path(const Path& path, const Position& pos);
// divide the path at a path-relative offset as measured in to_length from start
pair<Path, Path> cut_path(const Path& path, size_t offset);
bool maps_to_node(const Path& p, id_t id);
// the position that starts just after the path ends
Position path_start_position(const Path& path);
Position path_end_position(const Path& path);
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

/// Switches the node ids in the path to the ones indicated by the translator
void translate_node_ids(Path& path, const unordered_map<id_t, id_t>& translator);
/// Replaces the node IDs in the path with the ones indicated by the
/// translator. Supports a single cut node in the source graph, where the given
/// number of bases of the given node were removed from its left or right side
/// when making the source graph from the destination graph.
void translate_node_ids(Path& path, const unordered_map<id_t, id_t>& translator, id_t cut_node, size_t bases_removed, bool from_right);
/// Switches the node ids and orientations in the path to the ones indicated by the translator
void translate_oriented_node_ids(Path& path, const unordered_map<id_t, pair<id_t, bool>>& translator);
/// Switches node ids and orientations in the path to the ones indicated by the translator
void translate_oriented_node_ids(Path& path, const function<pair<id_t, bool>(id_t)>& translator);
    
// the first position on the path
pos_t initial_position(const Path& path);
// the last position on the path
pos_t final_position(const Path& path);

// Turn a list of node traversals into a path
Path path_from_node_traversals(const list<NodeTraversal>& traversals);

// Remove the paths with names matching the predicate from the graph.
// Store them in the list unless it is nullptr.
void remove_paths(Graph& graph, const function<bool(const string&)>& paths_to_take, std::list<Path>* matching);

// Get a Path from a handle graph
Path path_from_path_handle(const PathHandleGraph& graph, path_handle_t path_handle);

// Wrap a Path in an Alignment
Alignment alignment_from_path(const HandleGraph& graph, const Path& path);


/*
 * STL implementations of the protobuf object for use in in-memory operations
 */
class edit_t {
public:
    edit_t() noexcept : _from_length(0), _to_length(0), _sequence() {}
    edit_t(const edit_t& other) noexcept = default;
    edit_t(edit_t&& other) noexcept = default;
    ~edit_t() noexcept = default;
    edit_t& operator=(const edit_t& other) = default;
    edit_t& operator=(edit_t&& other) = default;
    inline int32_t from_length() const;
    inline void set_from_length(int32_t l);
    inline int32_t to_length() const;
    inline void set_to_length(int32_t l);
    inline const string& sequence() const;
    inline void set_sequence(const string& s);
    inline string* mutable_sequence();
    inline bool operator==(const edit_t& other) const;
    inline bool operator!=(const edit_t& other) const;
private:
    int32_t _from_length;
    int32_t _to_length;
    string _sequence;
};

// the mapping_t name is already taken
class path_mapping_t {
public:
    path_mapping_t() noexcept = default;
    path_mapping_t(const path_mapping_t& other) noexcept = default;
    path_mapping_t(path_mapping_t&& other) noexcept = default;
    ~path_mapping_t() noexcept = default;
    path_mapping_t& operator=(const path_mapping_t& other) = default;
    path_mapping_t& operator=(path_mapping_t&& other) = default;
    inline const position_t& position() const;
    inline position_t* mutable_position();
    inline const vector<edit_t>& edit() const;
    inline const edit_t& edit(size_t i) const;
    inline vector<edit_t>* mutable_edit();
    inline edit_t* mutable_edit(size_t i);
    inline edit_t* add_edit();
    inline size_t edit_size() const;
    inline bool operator==(const path_mapping_t& other) const;
    inline bool operator!=(const path_mapping_t& other) const;
private:
    position_t _position;
    vector<edit_t> _edit;
};

class path_t {
public:
    path_t() noexcept = default;
    path_t(const path_t& other) noexcept = default;
    path_t(path_t&& other) noexcept = default;
    ~path_t() noexcept = default;
    path_t& operator=(const path_t& other) = default;
    path_t& operator=(path_t&& other) = default;
    inline const vector<path_mapping_t>& mapping() const;
    inline const path_mapping_t& mapping(size_t i) const;
    inline vector<path_mapping_t>* mutable_mapping();
    inline path_mapping_t* mutable_mapping(size_t i);
    inline path_mapping_t* add_mapping();
    inline void clear_mapping();
    inline size_t mapping_size() const;
    inline bool operator==(const path_t& other) const;
    inline bool operator!=(const path_t& other) const;
private:
    vector<path_mapping_t> _mapping;
};

void from_proto_edit(const Edit& proto_edit, edit_t& edit);
void to_proto_edit(const edit_t& edit, Edit& proto_edit);
void from_proto_mapping(const Mapping& proto_mapping, path_mapping_t& mapping);
void to_proto_mapping(const path_mapping_t& mapping, Mapping& proto_mapping);
void from_proto_path(const Path& proto_path, path_t& path);
void to_proto_path(const path_t& path, Path& proto_path);


// repeated functions for the new path_t class
void translate_node_ids(path_t& path, const unordered_map<id_t, id_t>& translator);
void translate_oriented_node_ids(path_t& path, const unordered_map<id_t, pair<id_t, bool>>& translator);
void translate_oriented_node_ids(path_t& path, const function<pair<id_t, bool>(id_t)>& translator);

int mapping_from_length(const path_mapping_t& mapping);
int path_from_length(const path_t& path);
int mapping_to_length(const path_mapping_t& mapping);
int path_to_length(const path_t& path);

path_mapping_t reverse_complement_mapping(const path_mapping_t& m,
                                          const function<int64_t(id_t)>& node_length);
path_t reverse_complement_path(const path_t& path,
                               const function<int64_t(id_t)>& node_length);
void reverse_complement_mapping_in_place(path_mapping_t* m,
                                         const function<int64_t(id_t)>& node_length);
void reverse_complement_path_in_place(path_t* path,
                                      const function<int64_t(id_t)>& node_length);

// the first position on the path
pos_t initial_position(const path_t& path);
// the last position on the path
pos_t final_position(const path_t& path);

int corresponding_to_length(const path_t& path, int from_length, bool from_end);
int corresponding_from_length(const path_t& path, int to_length, bool from_end);

string debug_string(const path_t& path);
string debug_string(const path_mapping_t& mapping);
string debug_string(const edit_t& edit);

/*
 * Implementations of inline methods
 */

/*
 * edit_t
 */
inline int32_t edit_t::from_length() const {
    return _from_length;
}
inline void edit_t::set_from_length(int32_t l) {
    _from_length = l;
}
inline int32_t edit_t::to_length() const {
    return _to_length;
}
inline void edit_t::set_to_length(int32_t l) {
    _to_length = l;
}
inline const string& edit_t::sequence() const {
    return _sequence;
}
inline void edit_t::set_sequence(const string& s) {
    _sequence = s;
}
inline string* edit_t::mutable_sequence() {
    return &_sequence;
}
inline bool edit_t::operator==(const edit_t& other) const {
    return (_to_length == other._to_length
            && _from_length == other._from_length
            && _sequence == other._sequence);
}
inline bool edit_t::operator!=(const edit_t& other) const {
    return !(*this == other);
}

/*
 * path_mapping_t
 */
inline const position_t& path_mapping_t::position() const {
    return _position;
}
inline position_t* path_mapping_t::mutable_position() {
    return &_position;
}
inline const vector<edit_t>& path_mapping_t::edit() const {
    return _edit;
}
inline const edit_t& path_mapping_t::edit(size_t i) const {
    return _edit[i];
}
inline vector<edit_t>* path_mapping_t::mutable_edit() {
    return &_edit;
}
inline edit_t* path_mapping_t::add_edit() {
    _edit.emplace_back();
    return &_edit.back();
}
inline edit_t* path_mapping_t::mutable_edit(size_t i) {
    return &_edit[i];
}
inline size_t path_mapping_t::edit_size() const {
    return _edit.size();
}
inline bool path_mapping_t::operator==(const path_mapping_t& other) const {
    return (_position == other._position && _edit == other._edit);
}
inline bool path_mapping_t::operator!=(const path_mapping_t& other) const {
    return !(*this == other);
}

/*
 * path_t
 */
inline const vector<path_mapping_t>& path_t::mapping() const {
    return _mapping;
}
inline const path_mapping_t& path_t::mapping(size_t i) const {
    return _mapping[i];
}
inline vector<path_mapping_t>* path_t::mutable_mapping() {
    return &_mapping;
}
inline path_mapping_t* path_t::mutable_mapping(size_t i) {
    return &_mapping[i];
}
inline path_mapping_t* path_t::add_mapping() {
    _mapping.emplace_back();
    return &_mapping.back();
}
inline void path_t::clear_mapping() {
    _mapping.clear();
}
inline size_t path_t::mapping_size() const {
    return _mapping.size();
}
inline bool path_t::operator==(const path_t& other) const {
    return _mapping == other._mapping;
}
inline bool path_t::operator!=(const path_t& other) const {
    return !(*this == other);
}
}

#endif

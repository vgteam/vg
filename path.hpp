#ifndef PATH_H
#define PATH_H

#include <iostream>
#include <algorithm>
#include <functional>
#include <set>
#include <list>
#include "json2pb.h"
#include "vg.pb.h"
#include "edit.hpp"
#include "hash_map.hpp"

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
    // for its path. The list in question is stroed above in _paths. Recall that
    // std::list iterators are bidirectional.
    map<Mapping*, list<Mapping>::iterator > mapping_itr;
    // This maps from Mapping* pointer to the name of the path it belongs to
    // (which can then be used to get the list its iterator belongs to).
    map<Mapping*, string> mapping_path;
    void sort_by_mapping_rank(void);
    void rebuild_mapping_aux(void);
    // ...we need this in order to get subsets of the paths in correct order
    map<Mapping*, size_t> mapping_path_order;
    // This maps from node ID to a set of path name, mapping instance pairs.
    // Note that we don't have a map for each node, because each node can appear
    // along any given path multiple times, with multiple Mapping* pointers.
    map<int64_t, map<string, set<Mapping*>>> node_mapping;
    
    void rebuild_node_mapping(void);
    //void sync_paths_with_mapping_lists(void);
    list<Mapping>::iterator remove_mapping(Mapping* m);
    list<Mapping>::iterator insert_mapping(list<Mapping>::iterator w,
                                           const string& path_name, const Mapping& m);
    void remove_paths(const set<string>& names);
    void keep_paths(const set<string>& name);
    void remove_node(int64_t id);
    bool has_path(const string& name);
    void to_json(ostream& out);
    list<Mapping>& get_path(const string& name);
    list<Mapping>& get_create_path(const string& name);
    list<Mapping>& create_path(const string& name);
    bool has_mapping(const string& name, const Mapping& m);
    bool has_node_mapping(int64_t id);
    bool has_node_mapping(Node* n);
    map<string, set<Mapping*>>& get_node_mapping(Node* n);
    map<string, set<Mapping*>>& get_node_mapping(int64_t id);
    // Go left along the path that this Mapping* belongs to, and return the
    // Mapping* there, or null if this Mapping* is the first in its path.
    Mapping* traverse_left(Mapping* mapping);
    // Go right along the path that this Mapping* belongs to, and return the
    // Mapping* there, or null if this Mapping* is the last in its path.
    Mapping* traverse_right(Mapping* mapping);
    // TODO: should this be a reference?
    string mapping_path_name(Mapping* m);
    set<string> of_node(int64_t id);
    bool are_consecutive_nodes_in_path(int64_t id1, int64_t id2, const string& path_name);
    size_t size(void) const;
    bool empty(void) const;
    void clear(void);
    void clear_node_ranks(void);
    //void add_node_mapping(Node* n);
    void load(istream& in);
    void write(ostream& out);
    void to_graph(Graph& g);
    // add mappings, assume sorted by default
    void append_mapping(const string& name, const Mapping& m);
    void append_mapping(const string& name, int64_t id, size_t rank = 0, bool is_reverse = false);
    void append(Paths& p);
    void append(Graph& g);
    void extend(Paths& p);
    void extend(const Path& p);
    void for_each(const function<void(Path&)>& lambda);
    void for_each_stream(istream& in, const function<void(Path&)>& lambda);
    void increment_node_ids(int64_t inc);
    // Replace the node IDs used as keys with those used as values.
    // This is only efficient to do in a batch.
    void swap_node_ids(hash_map<int64_t, int64_t> id_mapping);
    void for_each_mapping(const function<void(Mapping*)>& lambda);
};

Path& increment_node_mapping_ids(Path& p, int64_t inc);
Path& append_path(Path& a, const Path& b);
const Paths paths_from_graph(Graph& g);
void parse_region(const string& target, string& name, int64_t& start, int64_t& end);
int path_to_length(const Path& path);
int path_from_length(const Path& path);
int mapping_to_length(const Mapping& m);
int mapping_from_length(const Mapping& m);
Path merge_paths(const Path& path1, const Path& path2, int& kept_path1, int& kept_path2);
Position first_path_position(const Path& path);
Position last_path_position(const Path& path);
int to_length(const Mapping& m);
int from_length(const Mapping& m);
bool mapping_ends_in_deletion(const Mapping& m);
bool mapping_starts_in_deletion(const Mapping& m);
bool mapping_is_total_deletion(const Mapping& m);
// Reverse-complement a Mapping and all the Edits in it. A function to get node
// lengths is needed, because the mapping will need to count its position from
// the other end of the node.
Mapping reverse_mapping(const Mapping& m, function<int64_t(int64_t)>& node_length);
// Reverse-complement a Path and all the Mappings in it. A function to get node
// lengths is needed, because the mappings will need to count their positions
// from the other ends of their nodes.
Path reverse_path(const Path& path, function<int64_t(int64_t)>& node_length);
// Simplify the path for addition as new material in the graph. Remove any
// mappings that are merely single deletions, merge adjacent edits of the same
// type, strip leading and trailing deletion edits on mappings, and make sure no
// mappings have missing positions.
Path simplify(const Path& p);
Mapping simplify(const Mapping& m);
Mapping merge(const Mapping& m, const Mapping& n);
Path concat_paths(const Path& path1, const Path& path2);
// divide mapping at reference-relative position
pair<Mapping, Mapping> cut_mapping(const Mapping& m, const Position& pos);
// divide mapping at target-relative offset (as measured in to_length)
pair<Mapping, Mapping> cut_mapping(const Mapping& m, size_t offset);
// divide path at reference-relative position
pair<Path, Path> cut_path(const Path& path, const Position& pos);
// divide the path at a path-relative offset as measured in to_length from start
pair<Path, Path> cut_path(const Path& path, size_t offset);
bool maps_to_node(const Path& p, int64_t id);
// Find all the points at which a Path enters or leaves nodes in the graph. Adds
// them to the given map by node ID of sets of bases in the node that will need
// to become the starts of new nodes. Note that some breakpoints may be at 0 or
// the past-the-end position of the node in question, since we don't look at the
// graph here.
void find_breakpoints(const Path& path, map<int64_t, set<int64_t>>& breakpoints);
// the position that starts just after the path ends
Position path_start(const Path& path);
Position path_end(const Path& path);
bool adjacent_mappings(const Mapping& m1, const Mapping& m2);


}

#endif

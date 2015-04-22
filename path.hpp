#ifndef PATH_H
#define PATH_H

#include <iostream>
#include <algorithm>
#include <functional>
#include <set>
#include <list>
#include "pb2json.h"
#include "vg.pb.h"
#include "json.hpp"

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

    //::google::protobuf::RepeatedPtrField< ::vg::Path >* _paths;
    map<string, list<Mapping> > _paths;
    map<Mapping*, list<Mapping>::iterator> _mapping_itr;
    map<int64_t, set<pair<string, Mapping*> > > node_mapping;
    void rebuild_node_mapping(void);
    void sync_paths_with_mapping_lists(void);
    void remove_paths(const set<string>& names);
    void keep_paths(const set<string>& name);
    void remove_node(int64_t id);
    bool has_path(const string& name);
    list<Mapping>& get_path(const string& name);
    list<Mapping>& get_create_path(const string& name);
    list<Mapping>& create_path(const string& name);
    bool has_mapping(const string& name, const Mapping& m);
    set<pair<string, Mapping*> >& get_node_mapping(Node* n);
    set<pair<string, Mapping*> >& get_node_mapping(int64_t id);
    set<string> of_node(int64_t id);
    size_t size(void) const;
    void clear(void);
    //void add_node_mapping(Node* n);
    void load(istream& in);
    void write(ostream& out);
    void to_graph(Graph& g);
    void append_mapping(const string& name, const Mapping& m);
    void append_mapping(const string& name, int64_t id, bool is_reverse = false);
    void append(Paths& p);
    void append(Graph& g);
    void extend(Paths& p);
    void extend(Path& p);
    void for_each(function<void(Path&)>& lambda);
    void for_each_stream(istream& in, function<void(Path&)>& lambda);
    void increment_node_ids(int64_t inc);
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

}

#endif

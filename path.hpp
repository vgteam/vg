#ifndef PATH_H
#define PATH_H

#include <iostream>
#include <algorithm>
#include <functional>
#include <set>
#include "pb2json.h"
#include "vg.pb.h"
#include "json.hpp"

namespace vg {

using namespace std;

class Paths : public vector<Path> {
public:
    map<int64_t, Path*> path_by_id;
    map<string, Path*> path_by_name;
    map<int64_t, set<pair<Path*, Mapping*> > > node_mapping;
    bool has_path(int64_t id);
    bool has_path(const string& name);
    Path* get_path(int64_t id);
    Path* get_path(const string& name);
    Path* get_create_path(const string& name);
    set<pair<Path*, Mapping*> >& get_node_mapping(Node* n);
    set<pair<Path*, Mapping*> >& get_node_mapping(int64_t id);
    //void add_node_mapping(Node* n);
    int64_t max_path_id(void);
    void load(istream& in);
    void write(ostream& out);
    void from_graph(Graph& g);
    void to_graph(Graph& g);
    void append_mapping(Path& p, Mapping& m);
    void append(Paths& p);
    void append(Graph& g);
    void extend(Paths& p);
    void extend(Path& p);
    Path* create_path(const string& name);
    void for_each(function<void(Path&)>& lambda);
    void for_each_stream(istream& in, function<void(Path&)>& lambda);
    void increment_node_ids(int64_t inc);
    void increment_path_ids(int64_t inc);
    void assign_path_ids(void);
    void for_each_mapping(const function<void(Mapping*)>& lambda);
};

Path& increment_node_mapping_ids(Path& p, int64_t inc);
Path& append_path(Path& a, const Path& b);
const Paths paths_from_graph(Graph& g);

}

#endif

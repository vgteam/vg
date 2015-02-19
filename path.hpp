#ifndef PATH_H
#define PATH_H

#include <iostream>
#include <algorithm>
#include <functional>
#include "pb2json.h"
#include "vg.pb.h"
#include "json.hpp"

namespace vg {

using namespace std;

class Paths : public vector<Path> {
public:
    void load(istream& in);
    void write(ostream& out);
    void from_graph(Graph& g);
    void append(Paths& p);
    void append(Graph& g);
    void extend(const Paths& p);
    void extend(const Path& p);
    void for_each(function<void(Path&)>& lambda);
    void for_each_stream(istream& in, function<void(Path&)>& lambda);
    void increment_ids(int64_t inc);
    void for_each_mapping(const function<void(Mapping*)>& lambda);
};

Path& increment_node_mapping_ids(Path& p, int64_t inc);
Path& append_path(Path& a, Path& b);
const Paths paths_from_graph(Graph& g);

}

#endif

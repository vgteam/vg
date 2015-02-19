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
    void append(Paths& p);
    void extend(Paths& p);
    void for_each(function<void(Path&)>& lambda);
    void for_each_stream(istream& in, function<void(Path&)>& lambda);
    void increment_ids(int64_t inc);
    void for_each_mapping(const function<void(Mapping*)>& lambda);
};

Path& increment_node_mapping_ids(Path& p, int64_t inc);
Path& extend_path(Path& a, Path& b);

}

#endif

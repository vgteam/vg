#include "path.hpp"
#include "stream.hpp"


namespace vg {

void Paths::load(istream& in) {
    uint64_t count = 0;
    function<void(uint64_t)> handle_count = [this, &count](uint64_t c) { count = c; };
    function<void(Path&)> lambda = [this](Path& p) {
        this->push_back(p);
    };
    stream::for_each(in, lambda, handle_count);
}

void Paths::write(ostream& out) {
    function<Path(uint64_t)> lambda =
        [this](uint64_t i) -> Path {
        return this->at(i);
    };
    stream::write(out, size(), lambda);
}

void Paths::for_each(function<void(Path&)>& lambda) {
    std::for_each(begin(), end(), lambda);
}

void Paths::for_each_stream(istream& in, function<void(Path&)>& lambda) {
    uint64_t count = 0;
    function<void(uint64_t)> handle_count = [this, &count](uint64_t c) { count = c; };
    stream::for_each(in, lambda, handle_count);
}

void Paths::append(Paths& paths) {
    map<string, Path*> path_name;
    for (auto& path : *this) {
        path_name[path.name()] = &path;
    }
    for (auto& path : paths) {
        auto n = path_name.find(path.name());
        if (n != path_name.end()) {
            extend_path(*n->second, path);
        }
    }
}

/*
void Paths::to_json(ostream& out) {
    std::for_each(begin(), end(), [this, &out](Path& p) { to_json(out, p); });
}
*/

void Paths::increment_ids(int64_t inc) {
    std::for_each(begin(), end(), [inc](Path& p) {
            increment_node_mapping_ids(p, inc);
        });
}

void Paths::extend(Paths& p) {
    reserve(size() + p.size());
    insert(end(), p.begin(), p.end());
}

Path& increment_node_mapping_ids(Path& p, int64_t inc) {
    for (int i = 0; i < p.mapping_size(); ++i) {
        Mapping* mapping = p.mutable_mapping(i);
        mapping->set_node_id(mapping->node_id()+inc);
    }
    return p;
}

Path& extend_path(Path& a, Path& b) {
    a.mutable_mapping()->MergeFrom(b.mapping());
    return a;
}



}

#include "path.hpp"
#include "stream.hpp"


namespace vg {

void Paths::load(istream& in) {
    uint64_t count = 0;
    function<void(uint64_t)> handle_count = [this, &count](uint64_t c) { count = c; };
    function<void(Path&)> lambda = [this](Path& p) {
        this->extend(p);
    };
    stream::for_each(in, lambda, handle_count);
}

void Paths::write(ostream& out) {
    function<Path(uint64_t)> lambda =
        [this](uint64_t i) -> Path {
        return this->_paths->Get(i);
    };
    stream::write(out, _paths->size(), lambda);
}

void Paths::for_each(function<void(Path&)>& lambda) {
    std::for_each(_paths->begin(), _paths->end(), lambda);
}

void Paths::for_each_mapping(const function<void(Mapping*)>& lambda) {
    for (auto& path : *_paths) {
        for (int64_t i = 0; i < path.mapping_size(); ++i) {
            Mapping* m = path.mutable_mapping(i);
            lambda(m);
        }
    }
}

void Paths::for_each_stream(istream& in, function<void(Path&)>& lambda) {
    uint64_t count = 0;
    function<void(uint64_t)> handle_count = [this, &count](uint64_t c) { count = c; };
    stream::for_each(in, lambda, handle_count);
}

void Paths::from_graph(Graph& g) {
    for (int i = 0; i < g.path_size(); ++i) {
        Path* p = g.mutable_path(i);
        extend(*p);
    }
}

void Paths::to_graph(Graph& g) {
    std::for_each(_paths->begin(), _paths->end(), [this, &g](Path& p) {
            Path* np = g.add_path();
            *np = p;
        });
}

void Paths::extend(Path& p) {
    Path* path = get_create_path(p.name());
    *path = p;
    for (int i = 0; i < path->mapping_size(); ++i) {
        //map<int64_t, set<pair<Path*, Mapping*> > > node_mapping;
        Mapping* m = path->mutable_mapping(i);
        get_node_mapping(m->node_id()).insert(make_pair(path, m));
    }
}

Path* Paths::create_path(const string& name) {
    Path* path = _paths->Add();
    path->set_name(name);
    path_by_name[name] = path;
    return path;
}

void Paths::extend(Paths& p) {
    for (auto& path : *p._paths) {
        extend(path);
    }
}

void Paths::append(Graph& g) {
    for (int i = 0; i < g.path_size(); ++i) {
        Path* path = g.mutable_path(i);
        Path* p = get_create_path(path->name());
        append_path_cache_nodes(*p, *path);
    }
}

void Paths::append(Paths& paths) {
    for (auto& path : *paths._paths) {
        Path* p = get_create_path(path.name());
        append_path_cache_nodes(*p, path);
    }
}

void Paths::append_path_cache_nodes(Path& a, Path& b) {
    // iterate over all the added elements in path
    // add them to the node mapping
    for (int i = 0; i < b.mapping_size(); ++i) {
        Mapping* m = a.add_mapping();
        *m = *b.mutable_mapping(i);
        get_node_mapping(m->node_id()).insert(make_pair(&a, m));
    }
}

// problem, won't allow us to keep multiple identical mapping to the same node,
// as will happen with looping paths
bool Paths::has_mapping(const string& name, Mapping& m) {
    auto& node_mapping = get_node_mapping(m.node_id());
    for (auto& p : node_mapping) {
        Path* path = p.first;
        if (path->name() == name) {
            return true;
        }
    }
    return false;
}

void Paths::append_mapping(const string& name, Mapping& m) {
    // TODO fix interface... this requires the path to be in this container
    Path* pt = get_create_path(name);
    assert(pt); // guard for now
// TODO check if we have the mapping already ?
    if (!has_mapping(name, m)) {
        Mapping* nm = pt->add_mapping();
        *nm = m; // and copy over the details
        // add it to the node mappings
        auto& ms = get_node_mapping(m.node_id());
        ms.insert(make_pair(pt, nm));
    }
}

void Paths::append_mapping(const string& name, int64_t id, bool is_reverse) {
    Mapping m;
    m.set_node_id(id);
    if (is_reverse) m.set_is_reverse(true);
    append_mapping(name, m);
}

bool Paths::has_path(const string& name) {
    return path_by_name.find(name) != path_by_name.end();
}

void Paths::increment_node_ids(int64_t inc) {
    std::for_each(_paths->begin(), _paths->end(), [inc](Path& p) {
            increment_node_mapping_ids(p, inc);
        });
    rebuild_node_mapping();
}

void Paths::rebuild_node_mapping(void) {
    // starts with paths and rebuilds the index
    //map<int64_t, set<pair<Path*, Mapping*> > > node_mapping;
    node_mapping.clear();
    for (auto& path : *_paths) {
        for (int i = 0; i < path.mapping_size(); ++i) {
            //map<int64_t, set<pair<Path*, Mapping*> > > node_mapping;
            Mapping* m = path.mutable_mapping(i);
            get_node_mapping(m->node_id()).insert(make_pair(&path, m));
        }
    }
}

void Paths::remove_node(int64_t id) {
    node_mapping.erase(id);
}

size_t Paths::size(void) const {
    return _paths->size();
}

Path* Paths::get_path(const string& name) {
    return path_by_name[name];
}

void Paths::remove_paths(const set<string>& names) {
    for (auto& nm : node_mapping) {
        auto mapping_copy = nm.second;
        auto& mapping = nm.second;
        for (auto& p : mapping_copy) {
            if (names.count(p.first->name())) {
                mapping.erase(p);
            }
        }
    }
    ::google::protobuf::RepeatedPtrField< ::vg::Path > new_paths;
    for (auto& path : *_paths) {
        if (!names.count(path.name())) {
            Path* new_path = new_paths.Add();
            *new_path = path;
        }
    }
    _paths->Clear();
    _paths->CopyFrom(new_paths);
}

void Paths::keep_paths(const set<string>& names) {
    for (auto& nm : node_mapping) {
        auto mapping_copy = nm.second;
        auto& mapping = nm.second;
        for (auto& p : mapping_copy) {
            if (!names.count(p.first->name())) {
                mapping.erase(p);
            }
        }
    }
    ::google::protobuf::RepeatedPtrField< ::vg::Path > new_paths;
    for (auto& path : *_paths) {
        if (names.count(path.name())) {
            Path* new_path = new_paths.Add();
            *new_path = path;
        }
    }
    _paths->Clear();
    _paths->CopyFrom(new_paths);
}

Path* Paths::get_create_path(const string& name) {
    Path* path = NULL;
    if (!has_path(name)) {
        path = create_path(name);
    } else {
        path = get_path(name);
    }
    return path;
}

set<pair<Path*, Mapping*> >& Paths::get_node_mapping(int64_t id) {
    return node_mapping[id];
}

set<pair<Path*, Mapping*> >& Paths::get_node_mapping(Node* n) {
    return node_mapping[n->id()];
}

set<string> Paths::of_node(int64_t id) {
    set<string> path_names;
    auto& node_mapping = get_node_mapping(id);
    for (auto& p : node_mapping) {
        path_names.insert(p.first->name());
    }
    return path_names;
}

Path& increment_node_mapping_ids(Path& p, int64_t inc) {
    for (int i = 0; i < p.mapping_size(); ++i) {
        Mapping* mapping = p.mutable_mapping(i);
        mapping->set_node_id(mapping->node_id()+inc);
    }
    return p;
}

Path& append_path(Path& a, const Path& b) {
    a.mutable_mapping()->MergeFrom(b.mapping());
    return a;
}

void parse_region(const string& target, string& name, int64_t& start, int64_t& end) {
    start = -1;
    end = -1;
    size_t foundFirstColon = target.find(":");
    // we only have a single string, use the whole sequence as the target
    if (foundFirstColon == string::npos) {
        name = target;
    } else {
        name = target.substr(0, foundFirstColon);
	    size_t foundRangeDash = target.find("-", foundFirstColon);
        if (foundRangeDash == string::npos) {
            start = atoi(target.substr(foundFirstColon + 1).c_str());
            end = start;
        } else {
            start = atoi(target.substr(foundFirstColon + 1, foundRangeDash - foundRangeDash - 1).c_str());
            end = atoi(target.substr(foundRangeDash + 1).c_str());
        }
    }

}

int path_to_length(const Path& path) {
    int l = 0;
    for (const auto& m : path.mapping()) {
        l += mapping_to_length(m);
    }
    return l;
}

int path_from_length(const Path& path) {
    int l = 0;
    for (const auto& m : path.mapping()) {
        l += mapping_from_length(m);
    }
    return l;
}

int mapping_to_length(const Mapping& m) {
    int l = 0;
    for (int i = 0; i < m.edit_size(); ++i) {
        const Edit& e = m.edit(i);
        if (e.has_from_length() && !e.has_to_length()) {
            l += e.from_length();
        } else {
            l += e.to_length();
        }
    }
    return l;
}

int mapping_from_length(const Mapping& m) {
    int l = 0;
    for (int i = 0; i < m.edit_size(); ++i) {
        const Edit& e = m.edit(i);
        l += e.from_length();
    }
    return l;

}


}

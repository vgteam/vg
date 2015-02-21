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
        return this->at(i);
    };
    stream::write(out, size(), lambda);
}

void Paths::for_each(function<void(Path&)>& lambda) {
    std::for_each(begin(), end(), lambda);
}

void Paths::for_each_mapping(const function<void(Mapping*)>& lambda) {
    for (auto& path : *this) {
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
    std::for_each(begin(), end(), [this, &g](Path& p) {
            Path* np = g.add_path();
            *np = p;
        });
}

void Paths::extend(Path& p) {
    Path* path = create_path(p.name());
    *path = p;
    for (int i = 0; i < path->mapping_size(); ++i) {
        //map<int64_t, set<pair<Path*, Mapping*> > > node_mapping;
        Mapping* m = path->mutable_mapping(i);
        get_node_mapping(m->node_id()).insert(make_pair(path, m));
    }
}

Path* Paths::create_path(const string& name) {
    emplace_back();
    Path* path = &(this->back());
    path->set_id(max_path_id()+1);
    path_by_id[path->id()] = path;
    path->set_name(name);
    path_by_name[name] = path;
    return path;
}

void Paths::extend(Paths& p) {
    reserve(size() + p.size());
    for (auto& path : p) {
        extend(path);
    }
}

void Paths::append(Graph& g) {
    Paths paths;
    for (int i = 0; i < g.path_size(); ++i) {
        const Path& path = g.path(i);
        Path* p = get_create_path(path.name());
        append_path(*p, path);
    }
    append(paths);
}

void Paths::append(Paths& paths) {
    for (auto& path : paths) {
        Path* p = get_create_path(path.name());
        append_path(*p, path);
    }
}

void Paths::append_mapping(Path& p, Mapping& m) {
    // TODO fix interface... this requires the path to be in this container
    Path* pt = get_path(p.id());
    assert(pt); // guard for now
    Mapping* nm = pt->add_mapping();
    *nm = m; // and copy over the details
    // add it to the node mappings
    auto& ms = get_node_mapping(m.node_id());
    ms.insert(make_pair(pt, nm));
}

bool Paths::has_path(int64_t id) {
    return path_by_id.find(id) != path_by_id.end();
}

bool Paths::has_path(const string& name) {
    return path_by_name.find(name) != path_by_name.end();
}

/*
void Paths::to_json(ostream& out) {
    std::for_each(begin(), end(), [this, &out](Path& p) { to_json(out, p); });
}
*/

void Paths::increment_node_ids(int64_t inc) {
    std::for_each(begin(), end(), [inc](Path& p) {
            increment_node_mapping_ids(p, inc);
        });
}

void Paths::increment_path_ids(int64_t inc) {
    std::for_each(begin(), end(), [inc](Path& p) {
            p.set_id(p.id()+inc);
        });
}

void Paths::assign_path_ids(void) {
    int id = 0;
    std::for_each(begin(), end(), [this, &id](Path& p) {
            //cerr << "assigning path id " << p.name() << " = " << id << endl;
            p.set_id(id++);
            path_by_id[p.id()] = &p;
        });
}

int64_t Paths::max_path_id(void) {
    if (path_by_id.size()) {
        return path_by_id.rend()->first;
    } else {
        return -1;
    }
}

Path* Paths::get_path(int64_t id) {
    return path_by_id[id];
}

Path* Paths::get_path(const string& name) {
    return path_by_name[name];
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

/*
void Paths::add_node_mapping(Node* n) {
    auto& m = get_node_mapping(n->id());
    for (int i = 0; i < n->path_id_size(); ++i) {
        m.insert(n->path_id(i));
    }
}
*/

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

}

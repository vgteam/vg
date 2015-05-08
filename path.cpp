#include "path.hpp"
#include "stream.hpp"


namespace vg {

Paths::Paths(void) {
    // noop
}

void Paths::load(istream& in) {
    uint64_t count = 0;
    function<void(uint64_t)> handle_count = [this, &count](uint64_t c) { count = c; };
    function<void(Path&)> lambda = [this](Path& p) {
        this->extend(p);
    };
    stream::for_each(in, lambda, handle_count);
}

void Paths::write(ostream& out) {
    vector<string> path_names;
    for (auto& p : _paths) {
        const string& name = p.first;
        path_names.push_back(name);
    }
    function<Path(uint64_t)> lambda =
        [this, &path_names](uint64_t i) -> Path {
        list<Mapping>& mappings = _paths[path_names.at(i)];
        Path path;
        for (auto& m : mappings) {
            Mapping* nm = path.add_mapping();
            *nm = m;
        }
        return path;
    };
    stream::write(out, _paths.size(), lambda);
}

void Paths::to_graph(Graph& g) {
    for (auto& p : _paths) {
        const string& name = p.first;
        list<Mapping>& mappings = p.second;
        Path* path = g.add_path();
        path->set_name(name);
        for (auto& m : mappings) {
            Mapping* nm = path->add_mapping();
            *nm = m;
        }
    }
}

void Paths::for_each(function<void(Path&)>& lambda) {
    for (auto& p : _paths) {
        const string& name = p.first;
        list<Mapping>& mappings = p.second;
        Path path;
        path.set_name(name);
        for (auto& m : mappings) {
            Mapping* nm = path.add_mapping();
            *nm = m;
        }
        lambda(path);
    }
}

void Paths::for_each_mapping(const function<void(Mapping*)>& lambda) {
    for (auto& p : _paths) {
        list<Mapping>& path = p.second;
        for (auto& m : path) {
            lambda(&m);
        }
    }
}

void Paths::for_each_stream(istream& in, function<void(Path&)>& lambda) {
    uint64_t count = 0;
    function<void(uint64_t)> handle_count = [this, &count](uint64_t c) { count = c; };
    stream::for_each(in, lambda, handle_count);
}

void Paths::extend(Path& p) {
    list<Mapping>& path = get_create_path(p.name());
    for (int i = 0; i < p.mapping_size(); ++i) {
        const Mapping& m = p.mapping(i);
        path.push_back(m);
        get_node_mapping(m.node_id()).insert(make_pair(p.name(), &path.back()));
    }
}

void Paths::extend(Paths& p) {
    for (auto& l : p._paths) {
        const string& name = l.first;
        list<Mapping>& path = l.second;
        for (auto& m : path) {
            append_mapping(name, m);         
        }
    }
}

void Paths::append(Paths& paths) {
    for (auto& p : paths._paths) {
        const string& name = p.first;
        const list<Mapping>& path = p.second;
        for (auto& m : path) {
            append_mapping(name, m);
        }
    }
}

void Paths::append(Graph& g) {
    for (int i = 0; i < g.path_size(); ++i) {
        const Path& p = g.path(i);
        for (int j = 0; j < p.mapping_size(); ++j) {
            const Mapping& m = p.mapping(j);
            append_mapping(p.name(), m);
        }
    }
}

Path& append_path(Path& a, const Path& b) {
    a.mutable_mapping()->MergeFrom(b.mapping());
    return a;
}

// problem, won't allow us to keep multiple identical mapping to the same node,
// as will happen with looping paths
bool Paths::has_mapping(const string& name, const Mapping& m) {
    auto& node_mapping = get_node_mapping(m.node_id());
    for (auto& p : node_mapping) {
        const string& path_name = p.first;
        if (path_name == name) {
            return true;
        }
    }
    return false;
}

void Paths::append_mapping(const string& name, const Mapping& m) {
    // TODO fix interface... this requires the path to be in this container
    list<Mapping>& pt = get_create_path(name);
// TODO check if we have the mapping already ?
    if (!has_mapping(name, m)) {
        pt.push_back(m);
        Mapping* mp = &pt.back();
        // add it to the node mappings
        auto& ms = get_node_mapping(m.node_id());
        ms.insert(make_pair(name, mp));
        list<Mapping>::iterator mi = pt.end(); --mi;
        mapping_itr[mp] = mi;
        mapping_path[mp] = name;
    }
}

void Paths::append_mapping(const string& name, int64_t id, bool is_reverse) {
    Mapping m;
    m.set_node_id(id);
    if (is_reverse) m.set_is_reverse(true);
    append_mapping(name, m);
}

bool Paths::has_path(const string& name) {
    return _paths.find(name) != _paths.end();
}

void Paths::increment_node_ids(int64_t inc) {
    for (auto& p : _paths) {
        const string& name = p.first;
        list<Mapping>& path = p.second;
        for (auto& m : path) {
            m.set_node_id(m.node_id()+inc);
        }
    }
    rebuild_node_mapping();
}

void Paths::rebuild_node_mapping(void) {
    // starts with paths and rebuilds the index
    node_mapping.clear();
    for (auto& p : _paths) {
        const string& path_name = p.first;
        list<Mapping>& path = p.second;
        for (auto& m : path) {
            get_node_mapping(m.node_id()).insert(make_pair(path_name, &m));
        }
    }
}

void Paths::rebuild_mapping_aux(void) {
//    map<Mapping*, list<Mapping>::iterator> mapping_itr
    mapping_itr.clear();
    mapping_path.clear();
    for (auto& p : _paths) {
        const string& path_name = p.first;
        list<Mapping>& path = p.second;
        for (list<Mapping>::iterator i = path.begin(); i != path.end(); ++i) {
            mapping_itr[&*i] = i;
            mapping_path[&*i] = path_name;
        }
    }
}

void Paths::remove_node(int64_t id) {
    node_mapping.erase(id);
}

list<Mapping>::iterator Paths::remove_mapping(Mapping* m) {
    if (mapping_path.find(m) == mapping_path.end()) cerr << "no mapping" << endl;
    const string& path_name = mapping_path[m];
    int64_t id = m->node_id();
    auto& x = _paths[path_name];
    list<Mapping>::iterator p = _paths[path_name].erase(mapping_itr[m]);
    if (has_node_mapping(id)) {
        auto& node_path_mapping = get_node_mapping(id);
        node_path_mapping.erase(make_pair(path_name, m));
        if (node_path_mapping.empty()) node_mapping.erase(id);
    }
    mapping_path.erase(m);
    mapping_itr.erase(m);
    return p;
}

list<Mapping>::iterator Paths::insert_mapping(list<Mapping>::iterator w, const string& path_name, const Mapping& m) {
    auto px = _paths.find(path_name);
    list<Mapping>& path = px->second;
    list<Mapping>::iterator p;
    if (path.empty()) {
        path.push_front(m);
        p = path.begin();
    } else {
        p = path.insert(w, m);
    }
    get_node_mapping(m.node_id()).insert(make_pair(path_name, &*p));
    mapping_path[&*p] = path_name;
    mapping_itr[&*p] = p;
    return p;
}

void Paths::to_json(ostream& out) {
    function<void(Path&)> lambda = [this, &out](Path& p) {
        char *json2 = pb2json(p);
        out << json2 <<endl;
        free(json2);
    };
    for_each(lambda);
}

size_t Paths::size(void) const {
    return _paths.size();
}

void Paths::clear(void) {
    _paths.clear();
    node_mapping.clear();
}

list<Mapping>& Paths::get_path(const string& name) {
    return _paths[name];
}

void Paths::remove_paths(const set<string>& names) {
    for (auto& name : names) {
        _paths.erase(name);
    }
    rebuild_node_mapping();
    rebuild_mapping_aux();
}

void Paths::keep_paths(const set<string>& names) {
    set<string> to_remove;
    for (auto& p : _paths) {
        if (!names.count(p.first)) {
            to_remove.insert(p.first);
        }
    }
    remove_paths(to_remove);
}

list<Mapping>& Paths::create_path(const string& name) {
    return _paths[name];
}

list<Mapping>& Paths::get_create_path(const string& name) {
    if (!has_path(name)) {
        return create_path(name);
    } else {
        return get_path(name);
    }
}

bool Paths::has_node_mapping(int64_t id) {
    return node_mapping.find(id) != node_mapping.end();
}

bool Paths::has_node_mapping(Node* n) {
    return node_mapping.find(n->id()) != node_mapping.end();
}

set<pair<string, Mapping*> >& Paths::get_node_mapping(int64_t id) {
    return node_mapping[id];
}

set<pair<string, Mapping*> >& Paths::get_node_mapping(Node* n) {
    return node_mapping[n->id()];
}

string Paths::mapping_path_name(Mapping* m) {
    auto n = mapping_path.find(m);
    if (n == mapping_path.end()) {
        return "";
    } else {
        return n->second;
    }
}

set<string> Paths::of_node(int64_t id) {
    set<string> path_names;
    auto& node_mapping = get_node_mapping(id);
    for (auto& p : node_mapping) {
        path_names.insert(p.first);
    }
    return path_names;
}

bool Paths::are_consecutive_nodes_in_path(int64_t id1, int64_t id2, const string& path_name) {
    if (of_node(id1).count(path_name) && of_node(id2).count(path_name)) {
        auto& p1 = get_node_mapping(id1);
        auto& p2 = get_node_mapping(id2);
        // is p1 directly before p2?
        list<Mapping>::iterator i1, i2;
        // note that this will get the first mapping in each path, not an arbitrary one
        // (we can have looping paths, so there could be several mappings per path)
        for (auto& nm : p1) {
            if (nm.first == path_name) i1 = mapping_itr[nm.second];
        }
        for (auto& nm : p2) {
            if (nm.first == path_name) i2 = mapping_itr[nm.second];
        }
        ++i1; // increment the first node's mapping iterator
        if (i1 == i2) return true;
    }
    return false;
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

void path_into_mappings(const Path& path, map<int64_t, vector<Mapping> >& mappings) {
    for (int i = 0; i < path.mapping_size(); ++i) {
        const Mapping& m = path.mapping(i);
        mappings[m.node_id()].push_back(m);
    }
}

// designed for merging longer paths with substantial overlap into one long path
Path merge_paths(vector<Path>& paths) {
    Path result;
    map<int64_t, map<int32_t, vector<Path*> > > paths_by_position;
    for (auto& path: paths) {
        if (path.mapping_size() > 0) {
            const Mapping& mapping = path.mapping(0);
            paths_by_position[mapping.node_id()][mapping.offset()].push_back(&path);
        }
    }
    for (auto& n : paths_by_position) {
        for (auto& o : n.second) {
            for (auto& p : o.second) {
                Path path = *p;
                // scan forward in this path until we're past the soft clips
                const Edit& first_edit = path.mapping(0).edit(0);
                bool has_softclip = (first_edit.from_length() == 0 && first_edit.to_length() > 0);
                int merge_from = has_softclip ? 1 : 0;
                int64_t merge_at_id = path.mapping(merge_from).node_id();
                // scan backwards in the result path until we get to the same node/offset in this path
                int cut_result_at = 0;
                for (int i = result.mapping_size()-1; i > 0; --i) {
                    if (result.mapping(i).node_id() == merge_at_id) {
                        cut_result_at = i; break;
                    }
                }
                if (cut_result_at == 0) {
                    cerr << "cannot merge paths, no overlap" << endl;
                    exit(1);
                }
                int to_remove = result.mapping_size() - cut_result_at;
                for (int i = 0; i < to_remove; ++i) {
                    result.mutable_mapping()->RemoveLast();
                }
                // cut off the softclipped mapping if there is one
                if (merge_from > 0) path.mutable_mapping()->erase(path.mutable_mapping()->begin());
                // append this path to the result
                result.mutable_mapping()->MergeFrom(path.mapping());
            }
        }
    }
    return result;
}


}

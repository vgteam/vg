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
        get_node_mapping(m.position().node_id()).insert(make_pair(p.name(), &path.back()));
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
    auto& node_mapping = get_node_mapping(m.position().node_id());
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
        auto& ms = get_node_mapping(m.position().node_id());
        ms.insert(make_pair(name, mp));
        list<Mapping>::iterator mi = pt.end(); --mi;
        mapping_itr[mp] = mi;
        mapping_path[mp] = name;
    }
}

void Paths::append_mapping(const string& name, int64_t id, bool is_reverse) {
    Mapping m;
    m.mutable_position()->set_node_id(id);
    m.set_is_reverse(is_reverse);
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
            m.mutable_position()->set_node_id(m.position().node_id()+inc);
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
            get_node_mapping(m.position().node_id()).insert(make_pair(path_name, &m));
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
    int64_t id = m->position().node_id();
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
    get_node_mapping(m.position().node_id()).insert(make_pair(path_name, &*p));
    mapping_path[&*p] = path_name;
    mapping_itr[&*p] = p;
    return p;
}

void Paths::to_json(ostream& out) {
    function<void(Path&)> lambda = [this, &out](Path& p) {
        out << pb2json(p) <<endl;
    };
    for_each(lambda);
}

size_t Paths::size(void) const {
    return _paths.size();
}

bool Paths::empty(void) const {
    return _paths.size() == 0;
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

Mapping* Paths::traverse_left(Mapping* mapping) {
    // Get the iterator for this Mapping*
    list<Mapping>::iterator place = mapping_itr.at(mapping);
    
    // Get the path name for this Mapping*
    string path_name = mapping_path_name(mapping);
    
    // Get the list that the iterator is in
    list<Mapping>& path_list = _paths.at(path_name);
    
    // If we're already the beginning, return null.
    if(place == path_list.begin()) {
        return nullptr;
    }
    
    // Else walk left and return the address of the stored Mapping. std::list
    // iterators are bidirectional, so we will be able to do it.
    place--;
    return &(*place);
}

Mapping* Paths::traverse_right(Mapping* mapping) {
    // Get the iterator for this Mapping*
    list<Mapping>::iterator place = mapping_itr.at(mapping);
    
    // Get the path name for this Mapping*
    string path_name = mapping_path_name(mapping);
    
    // Get the list that the iterator is in
    list<Mapping>& path_list = _paths.at(path_name);
    
    // Advance the iterator right.
    place++;
    
    // If we're at the end, return null
    if(place == path_list.end()) {
        return nullptr;
    }
    
    // Else return the address of the stored Mapping.
    return &(*place);
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
        l += e.to_length();
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
        mappings[m.position().node_id()].push_back(m);
    }
}

int softclip_start(const Mapping& mapping) {
    int from_length = 0;
    int to_length = 0;
    int i = 0;
    while (i < mapping.edit_size() && from_length == 0) {
        from_length += mapping.edit(i).from_length();
        if (from_length > 0) break;
        to_length += mapping.edit(i).to_length();
        ++i;
    }
    return to_length;
}

int softclip_end(const Mapping& mapping) {
    int from_length = 0;
    int to_length = 0;
    int i = mapping.edit_size()-1;
    while (i >= 0 && from_length == 0) {
        from_length += mapping.edit(i).from_length();
        if (from_length > 0) break;
        to_length += mapping.edit(i).to_length();
        --i;
    }
    return to_length;
}

// returns the first non-softclip position in the path
Position first_path_position(const Path& path) {
    // step through soft clips
    int i = 0;
    while (i < path.mapping_size()) {
        if (from_length(path.mapping(i))) break;
        ++i;
    }
    if (i == path.mapping_size()) { cerr << "[vg::Path] end of path without a from_length" << endl; exit(1); }
    const Mapping& mapping = path.mapping(i);
    // find the soft clip length here
    Position pos = mapping.position();
    pos.set_offset(pos.offset()+softclip_start(mapping));
    return pos;
}

Position last_path_position(const Path& path) {
    int i = path.mapping_size()-1;
    while (i >= 0) {
        if (from_length(path.mapping(i))) break;
        --i;
    }
    if (i < 0) { cerr << "[vg::Path] start of path without a from_length" << endl; exit(1); }
    const Mapping& mapping = path.mapping(i);
    // find the soft clip length here
    Position pos = mapping.position();
    pos.set_offset(pos.offset() + from_length(mapping));
    return pos;
}

int to_length(const Mapping& m) {
    int l = 0;
    for (int i = 0; i < m.edit_size(); ++i) {
        const Edit& e = m.edit(i);
        l += e.to_length();
    }
    return l;
}

int from_length(const Mapping& m) {
    int l = 0;
    for (int i = 0; i < m.edit_size(); ++i) {
        const Edit& e = m.edit(i);
        l += e.from_length();
    }
    return l;
}

// designed for merging longer paths with substantial overlap into one long path
Path merge_paths(const Path& path1, const Path& path2, int& kept_path1, int& kept_path2) {

    // how much of each path do we keep?
    // used when manipulating path sequences elsewhere
    kept_path1 = 0;
    kept_path2 = 0;

    // todo.... handle edge case
    // where we match positions that split a single node
    // probably we are going to get bit by this

    int64_t p1_start = path1.mapping(0).position().node_id();
    int64_t p1_end = path1.mapping(path1.mapping_size()-1).position().node_id();
    int64_t p2_start = path2.mapping(0).position().node_id();
    int64_t p2_end = path2.mapping(path2.mapping_size()-1).position().node_id();
    //cerr << p1_start << "-" << p1_end << " should overlap " << p2_start << "-" << p2_end << endl;

    // scan through the second until we match the tail of the first
    int p1 = path1.mapping_size()-1;
    int p2 = 0; // and do the same for the first mapping
    //int64_t p1_node, p2_node;
    const Mapping* p1mp = &path1.mapping(p1);
    const Mapping* p2mp = &path2.mapping(p2);
    while (p2 < path2.mapping_size()) {
        p2mp = &path2.mapping(p2);
        if (p1mp->position().node_id() == p2mp->position().node_id()) break;
        ++p2;
        //cerr << p1mp->position().node_id() << " --- " << p2mp->position().node_id() << endl;
    }

    const Mapping& p1m = *p1mp;
    const Mapping& p2m = *p2mp;

    //cerr << p1m.position().node_id() << " vs " << p2m.position().node_id() << endl;
    assert(p1mp->position().node_id() == p2mp->position().node_id());

    Path result;

    for (int i = 0; i < p1; ++i) {
        Mapping* m = result.add_mapping();
        *m = path1.mapping(i);
        kept_path1 += to_length(*m);
    }

    // we are now pointing at the same node
    // however, the alignment may not completely overlap the node
    // so we need to ensure that we have overlap
    //if (p1 == path1.mapping_size()-1 || p2 == 0) {
        // we could be partially mapping to a node
        // add up the total sequence length across the two
        //to_length(
    if (p2m.position().offset()) {
        if (from_length(p1m) < p2m.position().offset()) {
            cerr << "[vg::Path] cannot merge paths as their ends do not overlap" << endl
                 << pb2json(p1m) << endl << pb2json(p2m) <<endl;
            exit(1);
        }
        // we have overlap or match
        kept_path1 += p2m.position().offset();
        kept_path2 += to_length(p2m);
        // make a new mapping that combines the two
        Mapping* m = result.add_mapping();
        *m->mutable_position() = p1m.position();
        // now add the edits from the previous mapping to the next mapping
        for (int i = 0; i < p1m.edit_size(); ++i) {
            Edit* e = m->add_edit();
            *e = p1m.edit(i);
        }
        // and skip this length in the next mapping
        int to_skip = mapping_to_length(*m) - p2m.position().offset();
        for (int i = 0; i < p2m.edit_size(); ++i) {
            if (to_skip > p2m.edit(i).to_length()) {
            } else {
                // divide the edit
                Edit* e = m->add_edit();
                *e = p2m.edit(i);
                e->set_to_length(e->to_length() - to_skip);
                e->set_from_length(e->from_length() - to_skip);
                if (!e->sequence().empty()) {
                    e->set_sequence(e->sequence().substr(to_skip));
                }
            }
        }
        // offset is 0
        m->mutable_position()->set_offset(0);
    } else {
        kept_path2 += to_length(p2m);
        Mapping* m = result.add_mapping();
        *m = p2m;
    }

    if (p2+1 < path2.mapping_size()) {
        for (int i = p2+1; i < path2.mapping_size(); ++i) {
            Mapping* m = result.add_mapping();
            *m = path2.mapping(i);
            kept_path2 += to_length(*m);
        }
    }

    return result;
}

Path simplify_deletions(const Path& p) {
    Path s;
    for (int i = 0; i < p.mapping_size(); ++i) {
        auto& m = p.mapping(i);
        if (m.edit_size() == 1 && edit_is_deletion(m.edit(0))) {
        } else {
            *s.add_mapping() = m;
        }
    }
    return s;
}

bool mapping_ends_in_deletion(const Mapping& m){
    return edit_is_deletion(m.edit(m.edit_size()-1));
}

bool mapping_starts_in_deletion(const Mapping& m) {
    return edit_is_deletion(m.edit(0));
}

bool mapping_is_total_deletion(const Mapping& m) {
    return m.edit_size() == 1 && edit_is_deletion(m.edit(0));
}



}

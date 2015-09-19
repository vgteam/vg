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

void Paths::extend(const Path& p) {
    const string& name = p.name();
    list<Mapping>& path = get_create_path(name);
    cerr << "extending with " << name << endl;
    for (int i = 0; i < p.mapping_size(); ++i) {
        const Mapping& m = p.mapping(i);
        append_mapping(name, m);
    }
    cerr << "done with extension" << endl;
}

// one of these should go away
void Paths::extend(Paths& p) {
    cerr << "in extend paths paths" << endl;
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
    // get or create the path with this name
    list<Mapping>& pt = get_create_path(name);
    // now if we haven't already supplied a mapping
    // add it
    if (!has_mapping(name, m)) {
        pt.push_back(m);
        Mapping* mp = &pt.back();
        // add it to the node mappings
        auto& ms = get_node_mapping(m.position().node_id());
        ms.insert(make_pair(name, mp));
        // and record its position in this list
        list<Mapping>::iterator mi = pt.end(); --mi;
        mapping_itr[mp] = mi;
        mapping_path[mp] = name;
    }
}

void Paths::append_mapping(const string& name, int64_t id, bool is_reverse) {
    cerr << "in weird append mapping" << endl;
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

void Paths::swap_node_ids(hash_map<int64_t, int64_t> id_mapping) {
    for (auto& p : _paths) {
        const string& name = p.first;
        list<Mapping>& path = p.second;
        for (auto& m : path) {
            // Look up the replacement ID
            auto replacement = id_mapping.find(m.position().node_id());
        
            if(replacement != id_mapping.end()) {
                // If there is a replacement, use it.
                m.mutable_position()->set_node_id((*replacement).second);
            }
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
    mapping_path_order.clear();
    for (auto& p : _paths) {
        const string& path_name = p.first;
        list<Mapping>& path = p.second;
        size_t order_in_path = 0;
        for (list<Mapping>::iterator i = path.begin(); i != path.end(); ++i) {
            mapping_itr[&*i] = i;
            mapping_path[&*i] = path_name;
            mapping_path_order[&*i] = order_in_path++;
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
    mapping_path_order.erase(m);
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

// concatenates paths
Path concat_paths(const Path& path1, const Path& path2) {
    Path res = path1;
    //cerr << "-------------------- concat thing ------------------" << endl;
    //cerr << pb2json(path1) << endl << pb2json(path2) << endl;
    // tack on the edits from the last
    auto& path1_back = path1.mapping(path1.mapping_size()-1);
    auto& path2_front = path2.mapping(0);

    // check if we have to splice the last mapping together
    if (!path2_front.has_position() || !path1_back.has_position() ||
        path1_back.position().node_id() == path2_front.position().node_id()) {
        auto* mapping = res.mutable_mapping(res.mapping_size()-1);
        // adapt unmapped paths (which look like insertions here)
        if (!path2_front.has_position() && path1_back.has_position()) {
            *mapping->mutable_position() = path1_back.position();
        } else if (!path1_back.has_position() && path2_front.has_position()) {
            *mapping->mutable_position() = path2_front.position();
        }
        // merge the edits from the second onto the last mapping
        for (size_t i = 0; i < path2_front.edit_size(); ++i) {
            *mapping->add_edit() = path2_front.edit(i);
        }
    } else {
        // just tack it on, it's on the next node
        *res.add_mapping() = path2_front;
    }
    // simply add the rest of the mappings
    for (size_t i = 1; i < path2.mapping_size(); ++i) {
        *res.add_mapping() = path2.mapping(i);
    }
    //cerr << ">>>>" << endl;
    //cerr << pb2json(res) << endl;
    //cerr << "<<<<<<<<<<<<<<<<<<<<<<<<<<<< end " << endl;
    return res;
}

// merge paths that overlap at a suitable position
// TODO probably buggy, moved to concatenate paths above for banded alignment
// this remains a useful method as it is Alignment-object independent and could work on graph paths
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

    //assert(p1mp->position().node_id() == p2mp->position().node_id());
    bool is_split = (p1mp->position().node_id() != p2mp->position().node_id());
    //if (is_split) cerr << "we have a split read deletion or jump" << endl;

    Path result;

    for (int i = 0; i < p1; ++i) {
        Mapping* m = result.add_mapping();
        *m = path1.mapping(i);
        kept_path1 += to_length(*m);
    }

    // we are now pointing at the same node
    if (!is_split && p2m.position().offset()) {
        // the mappings are against the same node
        // make a new mapping that combines the two
        Mapping* m = result.add_mapping();
        *m->mutable_position() = p1m.position();
        // now add the edits from the previous mapping to the next mapping
        for (int i = 0; i < p1m.edit_size(); ++i) {
            Edit* e = m->add_edit();
            *e = p1m.edit(i);
            kept_path1 += e->to_length();
        }
        // local deletion on this node
        if ((from_length(p1m) < p2m.position().offset())) {
            // add a deletion
            Edit* e = m->add_edit();
            e->set_from_length(p2m.position().offset() - from_length(p1m));
        }
        // and skip this length in the next mapping
        // caution, this is based on graph, not alignment coordinates
        int to_skip = max((int)0, (int)mapping_from_length(*m) - (int)p2m.position().offset());
        //cerr << "to_skip = " << to_skip << endl;
        size_t skipped = 0;
        size_t j = 0;
        for ( ; j < p2m.edit_size(); ++j) {
            auto& f = p2m.edit(j);
            if (f.from_length() + skipped > to_skip) {
                break;
            }
            skipped += f.from_length();
        }
        // now we're pointing at the edit to divide
        {
            auto& f = p2m.edit(j++);
            size_t skip_here = to_skip - skipped;
            Edit* e = m->add_edit();
            cerr << "inner old " << pb2json(f) << endl;
            *e = f;
            e->set_to_length(e->to_length() - skip_here);
            e->set_from_length(e->from_length() - skip_here);
            if (!e->sequence().empty()) {
                e->set_sequence(e->sequence().substr(skip_here));
            }
            //cerr << "inner new " << pb2json(*e) << endl;
            kept_path2 += e->to_length();
        }
        // now let's add in the rest of the edits
        for (size_t i = 0; i < p2m.edit_size(); ++i) {
            auto& e = p2m.edit(i);
            *m->add_edit() = e;
            kept_path2 += e.to_length();
        }
        
        for (int i = 0; i < p2m.edit_size(); ++i) {
            if (to_skip > skipped + p2m.edit(i).to_length()) {
            } else {
                // divide the edit
                Edit* e = m->add_edit();
                cerr << "inner old " << pb2json(p2m.edit(i)) << endl;
                *e = p2m.edit(i);
                e->set_to_length(e->to_length() - to_skip);
                e->set_from_length(e->from_length() - to_skip);
                if (!e->sequence().empty()) {
                    e->set_sequence(e->sequence().substr(to_skip));
                }
                //cerr << "inner new " << pb2json(*e) << endl;
                kept_path2 += e->to_length();
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

Path simplify(const Path& p) {
    Path s;
    s.set_name(p.name());
    for (int i = 0; i < p.mapping_size(); ++i) {
        auto& m = p.mapping(i);
        // remove wholly-deleted mappings as these are redundant
        if (m.edit_size() == 1 && edit_is_deletion(m.edit(0))) continue;
        Mapping n; // our new mapping for this node
        // always take the old position
        *n.mutable_position() = m.position();
        Edit e = m.edit(0);
        for (int j = 1; j < m.edit_size(); ++j) {
            auto& f = m.edit(j);
            // if the edit types are the same, merge them
            if (edit_is_match(e) && edit_is_match(f)
                || edit_is_sub(e) && edit_is_sub(f)
                || edit_is_deletion(e) && edit_is_deletion(f)
                || edit_is_insertion(e) && edit_is_insertion(f)) {
                // will be 0 for insertions, and + for the rest
                e.set_from_length(e.from_length()+f.from_length());
                // will be 0 for deletions, and + for the rest
                e.set_to_length(e.to_length()+f.to_length());
                // will be empty for both or have sequence for both
                e.set_sequence(e.sequence() + f.sequence());
            } else {
                // mismatched types are just put on
                *n.add_edit() = e;
                e = f;
            }
        }
        // and keep the last edit
        *n.add_edit() = e;
        *s.add_mapping() = n;
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

// ref-relative
pair<Mapping, Mapping> cut_mapping(const Mapping& m, const Position& pos) {
    Mapping left, right;
    cerr << "cutting mapping " << pb2json(m) << " at pos " << pb2json(pos) << endl;
    // left always has the position of the input mapping
    *left.mutable_position() = m.position();
    if (m.position().node_id() != pos.node_id()) {
        left = m;
    } else { // we're on the node
        if (pos.offset() == m.position().offset()) {
            // we will get a 0-length left
            cerr << "should get a 0-length left" << endl;
            right = m;
        } else if (pos.offset() >= mapping_from_length(m)) {
            cerr << "should get a 0-length right" << endl;
            // or a 0-length right
            left = m;
        } else {
            cerr << "we need to cut the mapping" << endl;
            // we need to cut the mapping
            // find the cut point and build the two mappings
            size_t seen = 0;
            size_t j = 0;
            // loop over those before our position
            for ( ; j < m.edit_size() && seen < pos.offset(); ++j) {
                auto& e = m.edit(j);
                if (seen + e.from_length() == pos.offset()) {
                    // this would be the last edit before the target position
                    // so we just drop it onto the last mapping of p1
                    cerr << "at last edit before the target position" << endl;
                    *left.add_edit() = e;
                } else if (seen + e.from_length() > pos.offset()) {
                    cerr << "we need to divide the edit" << endl;
                    // we need to divide this edit
                    auto edits = cut_edit_at_from(e, seen + e.from_length() - pos.offset());
                    cerr << "got edits " << pb2json(edits.first) << " and " << pb2json(edits.second) << endl;
                    *left.add_edit() = edits.first;
                    *right.mutable_position() = pos;
                    *right.add_edit() = edits.second;
                }
                seen += e.from_length();
            }
            // now we add to the second path
            assert(seen >= pos.offset());
            for ( ; j < m.edit_size(); ++j) {
                *right.add_edit() = m.edit(j);
            }
        }
    }
    return make_pair(left, right);
}

// mapping-relative
pair<Mapping, Mapping> cut_mapping(const Mapping& m, size_t offset) {
    Mapping left, right;
    // left always has the position of the input mapping
    *left.mutable_position() = m.position();
    // nothing to cut
    if (offset == 0) {
        // we will get a 0-length left
        right = m;
    } else if (offset >= mapping_to_length(m)) {
        // or a 0-length right
        left = m;
    } else {
        // we need to cut the mapping
        // find the cut point and build the two mappings
        size_t seen = 0;
        size_t j = 0;
        // loop over those before our position
        for ( ; j < m.edit_size() && seen < offset; ++j) {
            auto& e = m.edit(j);
            cerr << "at edit " << pb2json(e) << endl;
            if (seen + e.to_length() > offset) {
                // we need to divide this edit
                auto edits = cut_edit_at_to(e, seen + e.to_length() - offset);
                *left.add_edit() = edits.first;
                if (m.has_position()) {
                    right.mutable_position()->set_node_id(m.position().node_id());
                    right.mutable_position()->set_offset(left.position().offset()
                                                         + mapping_from_length(left));
                }
                *right.add_edit() = edits.second;
            } else {
                // this would be the last edit before the target position
                // so we just drop it onto the last mapping of p1
                *left.add_edit() = e;
            }
            seen += e.to_length();
        }
        // now we add to the second path
        assert(seen >= offset);
        for ( ; j < m.edit_size(); ++j) {
            *right.add_edit() = m.edit(j);
        }
    }
    cerr << "cut mappings " << endl
         << "------left " << pb2json(left) << endl << "------right " << pb2json(right) << endl;
    return make_pair(left, right);
}

// divide path at reference-relative position
pair<Path, Path> cut_path(const Path& path, const Position& pos) {
    Path p1, p2;
    size_t i = 0;
    // seek forward to the cut point
    for ( ; i < path.mapping_size(); ++i) {
        auto& m = path.mapping(i);
        cerr << "seeking cut pos " << pb2json(pos) << " at mapping " << pb2json(m) << endl;
        // the position is in this node, so make the cut
        if (m.position().node_id() == pos.node_id()) {
            cerr << "making cuts" << endl;
            auto mappings = cut_mapping(m, pos);
            cerr << "left cut " << pb2json(mappings.first) << " and right " << pb2json(mappings.second) << endl;
            // and save the cuts
            *p1.add_mapping() = mappings.first;
            *p2.add_mapping() = mappings.second;
            ++i; // we don't increment our mapping index when we break here
            break;
        } else {
            // otherwise keep adding the mappings onto our first path
            *p1.add_mapping() = m;
        }
    }
    // add in the rest of the edits
    for ( ; i < path.mapping_size(); ++i) {
        auto& m = path.mapping(i);
        *p2.add_mapping() = m;
    }
    cerr << "---cut_path left " << pb2json(p1) << endl << "---and right " << pb2json(p2) << endl;
    return make_pair(p1, p2);
}

// divide the path at a path-relative offset as measured in to_length from start
pair<Path, Path> cut_path(const Path& path, size_t offset) {
    Path p1, p2;
    size_t seen = 0;
    size_t i = 0;
    // seek forward to the cut point
    for ( ; i < path.mapping_size() && seen < offset; ++i) {
        auto& m = path.mapping(i);
        // the position is in this node, so make the cut
        if (seen + mapping_to_length(m) == offset) {
            *p1.add_mapping() = m;
        } else if (seen + mapping_to_length(m) > offset) {
            auto mappings = cut_mapping(m, offset - seen);
            // and save the cuts
            *p1.add_mapping() = mappings.first;
            *p2.add_mapping() = mappings.second;
            ++i; // we don't increment our mapping index when we break here
            seen += mapping_to_length(m); // same problem
            break;
        } else {
            // otherwise keep adding the mappings onto our first path
            *p1.add_mapping() = m;
        }
        seen += mapping_to_length(m);
    }
    cerr << "seen " << seen << " offset " << offset << endl;
    assert(seen >= offset);
    // add in the rest of the edits
    for ( ; i < path.mapping_size(); ++i) {
        auto& m = path.mapping(i);
        *p2.add_mapping() = m;
    }
    cerr << "---cut_path left " << pb2json(p1) << endl << "---and right " << pb2json(p2) << endl;
    return make_pair(p1, p2);
}

bool maps_to_node(const Path& p, int64_t id) {
    for (size_t i = 0; i < p.mapping_size(); ++i) {
        if (p.mapping(i).position().node_id() == id) return true;
    }
    return false;
}

}

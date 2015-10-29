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
    for (int i = 0; i < p.mapping_size(); ++i) {
        const Mapping& m = p.mapping(i);
        append_mapping(name, m);
    }
}

// one of these should go away
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
    sort_by_mapping_rank();
    rebuild_mapping_aux();
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
    if (node_mapping.find(name) == node_mapping.end()) return false; // no mappings for path
    for (auto* mp : node_mapping[name]) {
        if (m.rank() == mp->rank()
            && m.is_reverse() == mp->is_reverse()) {
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
        ms[name].insert(mp);
        // and record its position in this list
        list<Mapping>::iterator mi = pt.end(); --mi;
        mapping_itr[mp] = mi;
        mapping_path[mp] = name;
        mapping_path_order[mp] = m.rank();
    }
}

void Paths::append_mapping(const string& name, int64_t id, size_t rank, bool is_reverse) {
    Mapping m;
    m.mutable_position()->set_node_id(id);
    m.set_is_reverse(is_reverse);
    if (rank) {
        m.set_rank(rank);
    } else {
        m.set_rank(get_path(name).size()+1); // rank is 1-based
    }
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
            get_node_mapping(m.position().node_id())[path_name].insert(&m);
        }
    }
}

// attempt to sort the paths based on the recorded ranks of the mappings
void Paths::sort_by_mapping_rank(void) {
    for (auto p = _paths.begin(); p != _paths.end(); ++p) {
        list<Mapping>& path = p->second;
        path.sort([](const Mapping& m1, const Mapping& m2) {
                return m1.rank() < m2.rank();
            });
    }
}

void Paths::rebuild_mapping_aux(void) {
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
            // if we have a rank already, use it
            if (i->rank()) {
                mapping_path_order[&*i] = i->rank();
            } else {
                // otherwise we set the rank based on what we've built
                mapping_path_order[&*i] = order_in_path;
                i->set_rank(order_in_path+1);
            }
            ++order_in_path;
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
        node_path_mapping[path_name].erase(m);
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
    get_node_mapping(m.position().node_id())[path_name].insert(&*p);
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

void Paths::clear_node_ranks(void) {
    for (auto p = _paths.begin(); p != _paths.end(); ++p) {
        list<Mapping>& path = p->second;
        for (auto m = path.begin(); m != path.end(); ++m) {
            Mapping& mapping = *m;
            mapping.clear_rank();
        }
    }
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

map<string, set<Mapping*>>& Paths::get_node_mapping(int64_t id) {
    return node_mapping[id];
}

map<string, set<Mapping*>>& Paths::get_node_mapping(Node* n) {
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
        vector<list<Mapping>::iterator> i1s, i2s;
        // note that this will get the first mapping in each path, not an arbitrary one
        // (we can have looping paths, so there could be several mappings per path)
        for (auto& mp : p1[path_name]) {
            i1s.push_back(mapping_itr[mp]);
        }
        for (auto& mp : p2[path_name]) {
            i2s.push_back(mapping_itr[mp]);
        }
        for (auto i1 : i1s) {
            ++i1; // increment the first node's mapping iterator
            for (auto i2 : i2s) {
                if (i1 == i2) return true;
            }
        }
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
    // move insertions from the front of the second path onto the back of the first
    // this does not change positions or anything else

    // We can merge two mappings to the same node only if they're in the same
    // direction and at compatible positions on the node. TODO: determine that
    // and merge mappings.

    // check if we have to splice the last mapping together
    if (!path2_front.has_position() || !path1_back.has_position()) {
        auto* mapping = res.mutable_mapping(res.mapping_size()-1);
        // adapt unmapped paths (which look like insertions here)
        if (!path2_front.has_position() && path1_back.has_position()) {
            // If one mapping has no position, it can't reference reference sequence
            assert(mapping_from_length(path2_front) == 0);
            *mapping->mutable_position() = path1_back.position();
            // Copy the reverse flag
            mapping->set_is_reverse(path1_back.is_reverse());
        } else if (!path1_back.has_position() && path2_front.has_position()) {
            // If one mapping has no position, it can't reference reference sequence
            assert(mapping_from_length(path1_back) == 0);
            *mapping->mutable_position() = path2_front.position();
            // Copy the reverse flag
            mapping->set_is_reverse(path2_front.is_reverse());
        }
        // merge the edits from the second onto the last mapping
        for (size_t i = 0; i < path2_front.edit_size(); ++i) {
            *mapping->add_edit() = path2_front.edit(i);
        }
    } else {
        // just tack it on, it's on the next node
        *res.add_mapping() = path2_front;
    }
    // add the rest of the mappings
    for (size_t i = 1; i < path2.mapping_size(); ++i) {
        *res.add_mapping() = path2.mapping(i);
    }
    if (path_from_length(res) != path_from_length(path1) + path_from_length(path2)
        || path_to_length(res) != path_to_length(path1) + path_to_length(path2)) {
        cerr << "error: concatenate fails to produce a path with from_length and to_length "
             << "equal to the sum of those of its inputs" << endl
             << "path1  " << pb2json(path1) << endl
             << "path2  " << pb2json(path1) << endl
             << "return " << pb2json(res) << endl;
        exit(1);
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
            //cerr << "inner old " << pb2json(f) << endl;
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
                //cerr << "inner old " << pb2json(p2m.edit(i)) << endl;
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
    for (size_t i = 0; i < p.mapping_size(); ++i) {
        auto m = p.mapping(i);
        
        // remove wholly-deleted mappings as these are redundant
        if (m.edit_size() == 1 && edit_is_deletion(m.edit(0))) continue;
        
        // if this isn't the first mapping
        // split off any insertions from the start
        // and push them to the last mapping
        if (i > 0) {
            size_t ins_at_start = 0;
            for (size_t j = 0; j < m.edit_size(); ++j) {
                auto& e = m.edit(j);
                if (!edit_is_insertion(e)) break;
                ins_at_start += e.to_length();
            }
            // if there are insertions at the start, move them left
            if (ins_at_start) {
                auto p = cut_mapping(m, ins_at_start);
                auto& ins = p.first;
                // cerr << "insertion " << pb2json(ins) << endl;
                // take the position from the original mapping
                m = p.second;
                *m.mutable_position() = ins.position();
                // cerr << "before and after " << pb2json(ins) << " and " << pb2json(m) << endl;
                Mapping* l = s.mutable_mapping(s.mapping_size()-1);
                for (size_t j = 0; j < ins.edit_size(); ++j) {
                    auto& e = ins.edit(j);
                    *l->add_edit() = e;
                }
            }
        }

        // handle the rest of path
        Mapping n_base; // our new mapping for this node
        Mapping* n = &n_base;
        // if we don't have a position, try to use the last mapping
        if (!m.has_position() || m.position().node_id()==0) {
            // We can't have a mapping referencing bases that don't exist.
            // Otherwise if we pull its stuff into the previous mapping, it can
            // run off the end of the node.
            assert(mapping_from_length(m) == 0);
            if (i > 0) {
                n = s.mutable_mapping(s.mapping_size()-1);
            } else {
                //cerr << "warning: path has no position in first mapping" << endl;
            }
        } else {
            // take the old position
            *n->mutable_position() = m.position();
            // Copy the is_reverse flag
            n->set_is_reverse(m.is_reverse());
        }
        
        size_t j = 0;
        // to simplify, we skip deletions
        // these are implied by jumps in the path
        for ( ; j < m.edit_size(); ++j) {
            if (!edit_is_deletion(m.edit(j))) {
                break;
            } else {
                // Adjust the offset by the size of the deletion. If we're going
                // forward on the node, this moves the mapping offset positive.
                // Otherwise it moves the mapping offset negative.
                n->mutable_position()->set_offset(n->position().offset() 
                    + m.edit(j).from_length() * (n->is_reverse() ? -1 : 1));
            }
        }
        if (j < m.edit_size()) {
            Edit e = m.edit(j++);
            for ( ; j < m.edit_size(); ++j) {
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
                    *n->add_edit() = e;
                    e = f;
                }
            }
            // and keep the last edit
            // if it isn't a deletion
            if (!edit_is_deletion(e)) *n->add_edit() = e;
        }
        // and store the mapping
        *s.add_mapping() = *n;
        
    }
    return s;    
}

bool mapping_ends_in_deletion(const Mapping& m){
    return m.edit_size() >= 1 && edit_is_deletion(m.edit(m.edit_size()-1));
}

bool mapping_starts_in_deletion(const Mapping& m) {
    return m.edit_size() >= 1 && edit_is_deletion(m.edit(0));
}

bool mapping_is_total_deletion(const Mapping& m) {
    return m.edit_size() == 1 && edit_is_deletion(m.edit(0));
}

Mapping reverse_mapping(const Mapping& m, function<int64_t(int64_t)>& node_length) {
    // Make a new reversed mapping
    Mapping reversed = m;
    
    if(m.has_position() && m.position().node_id() != 0) {
        // Since we're flipping m, which way do we move its start position in
        // the coordinates of its node? Towards the node's start if we were
        // reverse, and away from the node's start if we were forward.
        int64_t direction = m.is_reverse() ? -1 : 1;
    
        // We need to have an offset from the same end of the node to the
        // *other* end of the mapping. For this we need the mapping's
        // from_length and the direction we need to move the mapping start point
        // (towards or away from the start of the node). Importantly, we do not
        // need the actual length of the node, since all offsets are from the
        // node's start, in node-local coordinates.
        reversed.mutable_position()->set_offset(m.position().offset() + (mapping_from_length(m) - 1) * direction);
    }
    
    
    // Flip the flag
    reversed.set_is_reverse(!m.is_reverse());
    
    // Clear out all the edits. TODO: we wasted time copying them
    reversed.clear_edit();
    
    for(size_t i = m.edit_size() - 1; i != (size_t) -1; i--) {
        // For each edit in reverse order, put it in reverse complemented
        *reversed.add_edit() = reverse_edit(m.edit(i));
    }
    
    return reversed;
}

Path reverse_path(const Path& path, function<int64_t(int64_t)>& node_length) {
    // Make a new reversed path
    Path reversed = path;
    
    // Clear out all the mappings. TODO: we wasted time copying them
    reversed.clear_mapping();
    
    for(size_t i = path.mapping_size() - 1; i != (size_t) -1; i--) {
        // For each mapping in reverse order, put it in reverse complemented and
        // measured from the other end of the node.
        *reversed.add_mapping() = reverse_mapping(path.mapping(i), node_length);
    }
    
    return reversed;
}

// ref-relative
pair<Mapping, Mapping> cut_mapping(const Mapping& m, const Position& pos) {
    Mapping left, right;
    assert(m.has_position() && m.position().node_id());
    
    // TODO: support reverse mappings
    assert(!m.is_reverse());
    
    //cerr << "cutting mapping " << pb2json(m) << " at pos " << pb2json(pos) << endl;
    // left always has the position of the input mapping
    *left.mutable_position() = m.position();
    if (m.position().node_id() != pos.node_id()) {
        left = m;
    } else { // we're on the node
        if (pos.offset() == m.position().offset()) {
            // we will get a 0-length left
            //cerr << "should get a 0-length left" << endl;
            right = m;
        } else if (pos.offset() >= mapping_from_length(m)) {
            //cerr << "should get a 0-length right" << endl;
            // or a 0-length right
            left = m;
        } else {
            //cerr << "we need to cut the mapping" << endl;
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
                    //cerr << "at last edit before the target position" << endl;
                    *left.add_edit() = e;
                } else if (seen + e.from_length() > pos.offset()) {
                    //cerr << "we need to divide the edit" << endl;
                    // we need to divide this edit
                    auto edits = cut_edit_at_from(e, seen + e.from_length() - pos.offset());
                    //cerr << "got edits " << pb2json(edits.first) << " and " << pb2json(edits.second) << endl;
                    *left.add_edit() = edits.first;
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
    *right.mutable_position() = pos;
    assert(!m.has_position()
           || (left.has_position()
               && left.position().node_id()
               && right.has_position()
               && right.position().node_id()));
    return make_pair(left, right);
}

// mapping-relative
pair<Mapping, Mapping> cut_mapping(const Mapping& m, size_t offset) {
    Mapping left, right;
    
    // both result mappings will be in the same orientation as the original
    left.set_is_reverse(m.is_reverse());
    right.set_is_reverse(m.is_reverse());
    
    //assert(m.has_position() && m.position().node_id());
    // left always has the position of the input mapping
    if (m.has_position()) *left.mutable_position() = m.position();
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
            //cerr << "at edit " << pb2json(e) << endl;
            if (seen + e.to_length() > offset) {
                // we need to divide this edit
                auto edits = cut_edit_at_to(e, seen + e.to_length() - offset);
                *left.add_edit() = edits.first;
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
    if (m.has_position()) {
        // The right mapping has a position on this same node
        right.mutable_position()->set_node_id(m.position().node_id());
        // The position is closer to the node start if we're looking at the node
        // in reverse, and further from the node start if we're looking at it
        // normally. We offset by the left mapping's size.
        right.mutable_position()->set_offset(left.position().offset()
                                             + mapping_from_length(left) * (m.is_reverse() ? -1 : 1));
    }
    assert(!m.has_position()
           || (left.has_position()
               && left.position().node_id()
               && right.has_position()
               && right.position().node_id()));
    //cerr << "cut mappings " << endl
    //     << "------left " << pb2json(left) << endl << "------right " << pb2json(right) << endl;
    return make_pair(left, right);
}

// divide path at reference-relative position
pair<Path, Path> cut_path(const Path& path, const Position& pos) {
    Path p1, p2;
    size_t i = 0;
    // seek forward to the cut point
    for ( ; i < path.mapping_size(); ++i) {
        auto& m = path.mapping(i);
        //cerr << "seeking cut pos " << pb2json(pos) << " at mapping " << pb2json(m) << endl;
        // the position is in this node, so make the cut
        if (m.position().node_id() == pos.node_id()) {
            //cerr << "making cuts" << endl;
            auto mappings = cut_mapping(m, pos);
            //cerr << "left cut " << pb2json(mappings.first) << " and right " << pb2json(mappings.second) << endl;
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
    assert(p1.mapping(0).has_position()
           && p1.mapping(0).position().node_id()
           && p2.mapping(0).has_position()
           && p2.mapping(0).position().node_id());
    //cerr << "---cut_path left " << pb2json(p1) << endl << "---and right " << pb2json(p2) << endl;
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
#ifdef debug
        cerr << "seeking cut offset " << offset << " at mapping " << pb2json(m) << endl;
#endif
        // the position is in this node, so make the cut
        if (seen + mapping_to_length(m) == offset) {
            *p1.add_mapping() = m;
        } else if (seen + mapping_to_length(m) > offset) {
#ifdef debug
            cerr << "making cuts" << endl;
#endif
            auto mappings = cut_mapping(m, offset - seen);
            // and save the cuts
            *p1.add_mapping() = mappings.first;
            *p2.add_mapping() = mappings.second;
#ifdef debug
            cerr << "left cut " << pb2json(mappings.first) << " and right " << pb2json(mappings.second) << endl;
#endif
            ++i; // we don't increment our mapping index when we break here
            seen += mapping_to_length(m); // same problem
            break;
        } else {
            // otherwise keep adding the mappings onto our first path
            *p1.add_mapping() = m;
        }
        seen += mapping_to_length(m);
    }
#ifdef debug
    cerr << "seen " << seen << " offset " << offset << endl;
#endif
    assert(seen >= offset);
    // add in the rest of the edits
    for ( ; i < path.mapping_size(); ++i) {
        auto& m = path.mapping(i);
        *p2.add_mapping() = m;
    }
    assert(!path.mapping(0).has_position()
           || (p1.mapping(0).has_position()
               && p1.mapping(0).position().node_id()
               && p2.mapping(0).has_position()
               && p2.mapping(0).position().node_id()));
#ifdef debug
    cerr << "---cut_path left " << pb2json(p1) << endl << "---and right " << pb2json(p2) << endl;
#endif
    return make_pair(p1, p2);
}

bool maps_to_node(const Path& p, int64_t id) {
    for (size_t i = 0; i < p.mapping_size(); ++i) {
        if (p.mapping(i).position().node_id() == id) return true;
    }
    return false;
}

void find_breakpoints(const Path& path, map<int64_t, set<int64_t>>& breakpoints) {
    // We need to work out what offsets we will need to break each node at, if
    // we want to add in all the new material and edges in this path.
    
    // TODO: Move into graph or give access to graph, so we can avoid making extra breakpoints
    
#ifdef debug
    cerr << "Processing path..." << endl;
#endif
    
    for (size_t i = 0; i < path.mapping_size(); ++i) {
        // For each Mapping in the path
        const Mapping& m = path.mapping(i);
        
        // What node are we on?
        int64_t node_id = m.position().node_id();
        
        if(node_id == 0) {
            // Skip Mappings that aren't actually to nodes.
            continue;
        }
        
        // See where the next edit starts in the node. It is always included
        // (even when the edit runs backward), unless the edit has 0 length in
        // the reference.
        int64_t edit_first_position = m.position().offset();
        
        // And in what direction we are moving (+1 or -1)
        int64_t direction = m.is_reverse() ? -1 : 1;
        
#ifdef debug
    cerr << "Processing mapping " << i << " to node " << node_id << " in direction " << direction <<
        " from " << edit_first_position << endl;
#endif
        
        for(size_t j = 0; j < m.edit_size(); ++j) {
            // For each Edit in the mapping
            const Edit& e = m.edit(j);
            
            // We know where the mapping starts in its node. But where does it
            // end (inclusive)? Note that if the edit has 0 reference length,
            // this may not actually be included in the edit (and
            // edit_first_position will be further along than
            // edit_last_position).
            int64_t edit_last_position = edit_first_position + (e.from_length() - 1) * direction;
            
#ifdef debug
            cerr << "Edit on " << node_id << " from " << edit_first_position << " to " << edit_last_position << endl;
            cerr << pb2json(e) << endl;
#endif 
            
            if(!edit_is_match(e) || (j == 0 && i > 0)) {
                // If this edit is not a perfect match, or if this is the first
                // edit in this mapping and we had a previous mapping we may
                // need to connect to, we need to make sure we have a breakpoint
                // at the start of this edit.
                
#ifdef debug
                cerr << "Need to break " << node_id << " at edit lower end " <<
                    max(edit_first_position, edit_first_position - direction) << endl;
#endif
                
                // We need to snip between edit_first_position and edit_first_position - direction.
                // Note that it doesn't matter if we put breakpoints at 0 and 1-past-the-end; those will be ignored.
                breakpoints[node_id].insert(max(edit_first_position, edit_first_position - direction));
            }
            
            if(!edit_is_match(e) || (j == m.edit_size() - 1 && i < path.mapping_size() - 1)) {
                // If this edit is not a perfect match, or if it is the last
                // edit in a mapping and we have a subsequent mapping we might
                // need to connect to, make sure we have a breakpoint at the end
                // of this edit.
                
#ifdef debug
                cerr << "Need to break " << node_id << " at past edit upper end " <<
                    max(edit_last_position, edit_last_position + direction) << endl;
#endif
                
                // We also need to snip between edit_last_position and edit_last_position + direction.
                breakpoints[node_id].insert(max(edit_last_position, edit_last_position + direction));
            }
            
            // TODO: for an insertion or substitution, note that we need a new
            // node and two new edges.
            
            // TODO: for a deletion, note that we need an edge. TODO: Catch
            // and complain about some things we can't handle (like a path with
            // a leading/trailing deletion)? Or just skip deletions when wiring.
            
            // Use up the portion of the node taken by this mapping, so we know
            // where the next mapping will start.
            edit_first_position = edit_last_position + direction;
        }
    }
    
}

// returns the start position, or an empty position if the path has no mappings with positions
Position path_start(const Path& path) {
    for (size_t i = 0; i < path.mapping_size(); ++i) {
        auto& mapping = path.mapping(i);
        if (mapping.has_position()) return mapping.position();
    }
}

// determine the path end
Position path_end(const Path& path) {
    Position pos;
    if (!path.mapping_size()) return pos;
    auto& last = path.mapping(path.mapping_size()-1);
    pos = last.position();
    pos.set_offset(pos.offset()+mapping_from_length(last));
    return pos;
}

bool adjacent_mappings(const Mapping& m1, const Mapping& m2) {
    return abs(m1.rank() - m2.rank()) == 1;
}

}

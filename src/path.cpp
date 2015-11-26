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

void Paths::for_each(const function<void(Path&)>& lambda) {
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

void Paths::for_each_stream(istream& in, const function<void(Path&)>& lambda) {
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
    // re-sort?
    sort_by_mapping_rank();
    rebuild_mapping_aux();
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
    sort_by_mapping_rank();
    rebuild_mapping_aux();
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
            && m.position().is_reverse() == mp->position().is_reverse()) {
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
    m.mutable_position()->set_is_reverse(is_reverse);
    if (rank) {
        m.set_rank(rank);
    } else {
        m.set_rank(get_path(name).size()+1); // rank is 1-based
    }
    append_mapping(name, m);
}

// these will split a mapping into two
// NB: each submapping ends up with the same rank as the parent
// however, they will be ordered correctly in the path
// we will need to normalize path ranks to make this right
pair<Mapping*, Mapping*> Paths::divide_mapping(Mapping* m, const Position& pos) {
    // this is needed to split mappinsg during e.g. normalization
    // but still ensure that the mappings are out there
    // what do we do?
    // first we take the mapping and divide it as we do
    auto n = cut_mapping(*m, pos);
    return replace_mapping(m, n);
}

pair<Mapping*, Mapping*> Paths::divide_mapping(Mapping* m, size_t offset) {
    auto n = cut_mapping(*m, offset);
    return replace_mapping(m, n);
}

pair<Mapping*, Mapping*> Paths::replace_mapping(Mapping* m, pair<Mapping, Mapping> n) {
    // then we remove it from the node it's pointing to
    // and replace it with the other two mappings
    // we'll give them the same rank, but record them in the right order
    // this leaves an invalid graph
    // there are a few ways to fix this--- they involve changing the way we record ranks
    // but for now it's going to be simplest if the calling context manages this
    auto& path_name = mapping_path_name(m);
    n.first.set_rank(m->rank());
    n.second.set_rank(m->rank());
    // remove the mapping, getting an iterator pointing to the element that was after it
    auto i = remove_mapping(m);
    // the insertion will happen in reverse order
    // because insert puts the position before the iterator it's given
    // so we first insert the second element
    i = insert_mapping(i, path_name, n.second);
    // and then the second
    auto j = insert_mapping(i, path_name, n.second);
    // and we return them in proper order
    return make_pair(&*i, &*j);
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

// compact the ranks preserving the relative rank order
void Paths::compact_ranks(void) {
    // first ensure the storage order of the mappings is correct
    sort_by_mapping_rank();
    // clear the ranks
    clear_mapping_ranks();
    // and rebuild them and other aux data structures
    rebuild_node_mapping();
    rebuild_mapping_aux();
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
    mapping_path.clear();
    mapping_path_order.clear();
}

void Paths::clear_mapping_ranks(void) {
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

map<string, map<int, Mapping>> Paths::get_node_mapping_copies_by_rank(int64_t id) {
    map<string, map<int, Mapping>> by_ranks;
    for (auto& p : get_node_mapping(id)) {
        auto& name = p.first;
        auto& mp = p.second;
        for (auto* m : mp) by_ranks[name][m->rank()] = *m;
    }
    return by_ranks;
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

const string Paths::mapping_path_name(Mapping* m) {
    auto n = mapping_path.find(m);
    if (n == mapping_path.end()) {
        return "";
    } else {
        return n->second;
    }
}

map<string, int> Paths::of_node(int64_t id) {
    map<string, int> path_counts;
    auto& node_mapping = get_node_mapping(id);
    for (auto& p : node_mapping) {
        path_counts[p.first]++;
    }
    return path_counts;
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
            mapping->mutable_position()->set_is_reverse(path1_back.position().is_reverse());
        } else if (!path1_back.has_position() && path2_front.has_position()) {
            // If one mapping has no position, it can't reference reference sequence
            assert(mapping_from_length(path1_back) == 0);
            *mapping->mutable_position() = path2_front.position();
            // Copy the reverse flag
            mapping->mutable_position()->set_is_reverse(path2_front.position().is_reverse());
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
    //cerr << "simplifying " << pb2json(p) << endl;
    // loop over the mappings in the path, doing a few things
    // exclude mappings that are total deletions
    // when possible, merge a mapping with the previous mapping
    // push inserted sequences to the left
    for (size_t i = 0; i < p.mapping_size(); ++i) {
        auto m = simplify(p.mapping(i));
        // remove wholly-deleted mappings as these are redundant
        if (m.edit_size() == 1 && edit_is_deletion(m.edit(0))) continue;
        // if this isn't the first mapping
        if (i > 0) {
            // refer to the last mapping
            Mapping* l = s.mutable_mapping(s.mapping_size()-1);
            // split off any insertions from the start
            // and push them to the last mapping
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
                for (size_t j = 0; j < ins.edit_size(); ++j) {
                    auto& e = ins.edit(j);
                    *l->add_edit() = e;
                }
            }
            // if our last mapping has no position, but we do, merge
            if (!l->has_position() && m.has_position()) {
                *l->mutable_position() = m.position();
            // otherwise, if we end at exactly the start position of the next mapping, we can merge
            } else if (l->has_position() && m.has_position()
                       && l->position().node_id() == m.position().node_id()
                       && l->position().offset() + mapping_from_length(*l) == m.position().offset()) {
                // we can merge the current mapping onto the old one
                *l = merge_mappings(*l, m);
            } else {
                *s.add_mapping() = m;
            }
        } else {
            *s.add_mapping() = m;
        }
    }
    // now set ranks
    for (size_t i = 0; i < s.mapping_size(); ++i) {
        auto* m = s.mutable_mapping(i);
        m->set_rank(i+1);
    }
    //cerr << "simplified " << pb2json(s) << endl;
    return s;
}

// simple merge
Mapping merge_mappings(const Mapping& m, const Mapping& n) {
    Mapping c = m;
    // add the edits on
    for (size_t i = 0; i < n.edit_size(); ++i) {
        *c.add_edit() = n.edit(i);
    }
    // merge anything that's identical
    return simplify(c);
}

Mapping simplify(const Mapping& m) {
    Mapping n;
    if (m.rank()) n.set_rank(m.rank());
    //cerr << "pre simplify " << pb2json(m) << endl;
    // get the position
    if (!m.has_position() || m.position().node_id()==0) {
        // do nothing
    } else {
        // take the old position
        *n.mutable_position() = m.position();
    }

    size_t j = 0;
    // to simplify, we skip deletions at the very start of the node
    // these are implied by jumps in the path from other nodes
    if (m.has_position() && m.position().offset() == 0) {
        for ( ; j < m.edit_size(); ++j) {
            if (!edit_is_deletion(m.edit(j))) {
                break;
            } else {
                // Adjust the offset by the size of the deletion.
                n.mutable_position()->set_offset(n.position().offset()
                                                 + m.edit(j).from_length());
            }
        }
    }

    // now go through the rest of the edits and see if we can merge them
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
                *n.add_edit() = e;
                e = f;
            }
        }
        // and keep the last edit
        // if it isn't a deletion
        if (!edit_is_deletion(e)) *n.add_edit() = e;
    }
    //cerr << "post simplify " << pb2json(n) << endl;
    return n;
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

bool mapping_is_simple_match(const Mapping& m) {
    return m.edit_size() == 1 && edit_is_match(m.edit(0));
}

const string mapping_sequence(const Mapping& mp, const Node& n) {
    assert(mp.position().node_id() == n.id());
    auto& node_seq = n.sequence();
    string seq;
    // todo reverse the mapping
    function<int64_t(int64_t)> lambda = [&node_seq](int64_t i){return node_seq.size();};
    Mapping m = (mp.position().is_reverse() ? reverse_mapping(mp, lambda) : mp);
    // then edit in the forward direction (easier)
    // and, if the mapping is reversed, finally reverse-complement the result
    size_t t = 0;
    size_t f = 0;
    for (size_t i = 0; i < m.edit_size(); ++i) {
        auto& e = m.edit(i);
        if (edit_is_match(e)) {
            seq.append(node_seq.substr(f, e.from_length()));
        } else if (edit_is_sub(e)) {
            seq.append(e.sequence());
        } else if (edit_is_insertion(e)) {
            seq.append(e.sequence());
        } else if (edit_is_deletion(e)) {
            // no-op
        }
        t += e.to_length();
        f += e.to_length();
    }
    // perfect match
    if (m.edit_size() == 0) {
        seq = node_seq;
    }
    return (mp.position().is_reverse() ? reverse_complement(seq) : seq);
}

Mapping reverse_mapping(const Mapping& m, const function<int64_t(int64_t)>& node_length) {
    // Make a new reversed mapping
    Mapping reversed = m;

    // switching around to the reverse strand requires us to find the end of the mapping
    // on the forward and convert to the reverse --- or vice versa
    if(m.has_position() && m.position().node_id() != 0) {
        Position p = m.position();
        // set to the end of the mapping
        p.set_offset(p.offset() + mapping_from_length(m));
        // then flip the position onto the other side
        *reversed.mutable_position()
            = make_position(reverse(make_pos_t(p),
                                    node_length(m.position().node_id())));
    }
    
    // Clear out all the edits. TODO: we wasted time copying them
    reversed.clear_edit();
    
    for(size_t i = m.edit_size() - 1; i != (size_t) -1; i--) {
        // For each edit in reverse order, put it in reverse complemented
        *reversed.add_edit() = reverse_edit(m.edit(i));
    }
    
    return reversed;
}

Path reverse_path(const Path& path, const function<int64_t(int64_t)>& node_length) {
    // Make a new reversed path
    Path reversed = path;
    
    // Clear out all the mappings. TODO: we wasted time copying them
    reversed.clear_mapping();
    
    for(size_t i = path.mapping_size() - 1; i != (size_t) -1; i--) {
        // For each mapping in reverse order, put it in reverse complemented and
        // measured from the other end of the node.
        *reversed.add_mapping() = reverse_mapping(path.mapping(i), node_length);
    }
    for (size_t i = 0; i < path.mapping_size(); ++i) {
        reversed.mutable_mapping(i)->set_rank(i+1);
    }
    
    return reversed;
}

// ref-relative
pair<Mapping, Mapping> cut_mapping(const Mapping& m, const Position& pos) {
    Mapping left, right;
    assert(m.has_position() && m.position().node_id());
    
    // TODO: support reverse mappings
    assert(!m.position().is_reverse());
    
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
    left.mutable_position()->set_is_reverse(m.position().is_reverse());
    right.mutable_position()->set_is_reverse(m.position().is_reverse());
    
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
        right.mutable_position()->set_offset(left.position().offset()
                                             + mapping_from_length(left));
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

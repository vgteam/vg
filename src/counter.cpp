#include "counter.hpp"

namespace vg {

Counter::Counter(void) : xgidx(nullptr) { }

Counter::Counter(xg::XG* xidx) : xgidx(xidx) {
    coverage_dynamic = gcsa::CounterArray(xgidx->seq_length, 8);
}

Counter::~Counter(void) {
    close_edit_tmpfile();
    remove_edit_tmpfile();
}

void Counter::load_from_file(const string& file_name) {
    ifstream in(file_name);
    load(in);
    is_compacted = true;
}

void Counter::save_to_file(const string& file_name) {
    ofstream out(file_name);
    serialize(out);
}

void Counter::load(istream& in) {
    coverage_civ.load(in);
    edit_csa.load(in);
}

size_t Counter::serialize(std::ostream& out,
                          sdsl::structure_tree_node* s,
                          std::string name) {
    make_compact();
    sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(s, name, sdsl::util::class_name(*this));
    size_t written = 0;
    written += coverage_civ.serialize(out, child, "graph_coverage_" + name);
    written += edit_csa.serialize(out, child, "edit_csa_" + name);
    sdsl::structure_tree::add_size(child, written);
    return written;
}

void Counter::make_compact(void) {
    // pack the dynamic countarry and edit coverage into the compact data structure
    if (is_compacted) return;
    // sync edit file
    close_edit_tmpfile();
    // temporaries for construction
    size_t basis_length = coverage_dynamic.size();
    int_vector<> coverage_iv;
    util::assign(coverage_iv, int_vector<>(basis_length));
    for (size_t i = 0; i < coverage_dynamic.size(); ++i) {
        coverage_iv[i] = coverage_dynamic[i];
    }
    util::assign(coverage_civ, coverage_iv);
    construct(edit_csa, edit_tmpfile_name, 1);
    // construct the record marker bitvector
    remove_edit_tmpfile();
    is_compacted = true;
}

void Counter::make_dynamic(void) {
    if (!is_compacted) return;
    // unpack the compact represenation into the countarray
    assert(false); // not implemented
    is_compacted = false;
}

void Counter::ensure_edit_tmpfile_open(void) {
    if (!tmpfstream.is_open()) {
        string base = ".vg-counter";
        edit_tmpfile_name = tmpfilename(base);
        tmpfstream.open(edit_tmpfile_name);
    }
}

void Counter::close_edit_tmpfile(void) {
    if (tmpfstream.is_open()) {
        tmpfstream << delim; // pad
        tmpfstream.close();
    }
}

void Counter::remove_edit_tmpfile(void) {
    if (!edit_tmpfile_name.empty()) {
        std::remove(edit_tmpfile_name.c_str());
        edit_tmpfile_name.clear();
    }
}

void Counter::add(const Alignment& aln, bool record_edits) {
    // open tmpfile if needed
    ensure_edit_tmpfile_open();
    // count the nodes, edges, and edits
    for (auto& mapping : aln.path().mapping()) {
        if (!mapping.has_position()) continue;
        size_t i = position_in_basis(mapping.position());
        for (auto& edit : mapping.edit()) {
            if (edit_is_match(edit)) {
                if (mapping.position().is_reverse()) {
                    for (size_t j = 0; j < edit.from_length(); ++j) {
                        coverage_dynamic.increment(i-j);
                    }
                } else {
                    for (size_t j = 0; j < edit.from_length(); ++j) {
                        coverage_dynamic.increment(i+j);
                    }
                }
            } else if (record_edits) {
                // we represent things on the forward strand
                string pos_repr = pos_key(i);
                string edit_repr = edit_value(edit, mapping.position().is_reverse());
                tmpfstream << pos_repr << edit_repr;
            }
            if (mapping.position().is_reverse()) {
                i -= edit.from_length();
            } else {
                i += edit.from_length();
            }
        }
    }
}

// find the position on the forward strand in the sequence vector
size_t Counter::position_in_basis(const Position& pos) {
    // get position on the forward strand
    if (pos.is_reverse()) {
        return (int64_t)xg_node_start(pos.node_id(), xgidx)
            + (int64_t)reverse(pos, xg_node_length(pos.node_id(), xgidx)).offset() - 1;
    } else {
        return (int64_t)xg_node_start(pos.node_id(), xgidx) + (int64_t)pos.offset();
    }
}

string Counter::pos_key(size_t i) {
    Position pos;
    size_t offset = 2;
    pos.set_node_id(i+offset);
    string pos_repr;
    pos.SerializeToString(&pos_repr);
    return delim + escape_delim(pos_repr);
}

string Counter::edit_value(const Edit& edit, bool revcomp) {
    string edit_repr;
    if (revcomp) {
        reverse_complement_edit(edit).SerializeToString(&edit_repr);
    } else {
        edit.SerializeToString(&edit_repr);
    }
    return delim + escape_delim(edit_repr);
}

string Counter::escape_delim(const string& s) {
    string escaped; escaped.reserve(s.size());
    for (size_t i = 0; i < s.size(); ++i) {
        char c = s[i];
        escaped.push_back(c);
        if (c == delim) escaped.push_back(delim);
    }
    return escaped;
}

string Counter::unescape_delim(const string& s) {
    string unescaped; unescaped.reserve(s.size());
    for (size_t i = 0; i < s.size()-1; ++i) {
        char c = s[i];
        char d = s[i+1];
        if (c == delim && d == delim) {
            unescaped.push_back(delim);
        } else {
            unescaped.push_back(c);
            if (i == s.size()-2) unescaped.push_back(d);
        }
    }
    return unescaped;
}

vector<Edit> Counter::edits_at_position(size_t i) {
    vector<Edit> edits;
    if (i == 0) return edits;
    string key = pos_key(i);
    auto occs = locate(edit_csa, key);
    for (size_t i = 0; i < occs.size(); ++i) {
        //cout << occs[i] << endl;
        // walk from after the key to the next end-sep
        size_t b = occs[i] + key.size() + 1;
        size_t e = b;
        while (!(extract(edit_csa, e, e)[0] == delim && extract(edit_csa, e+1, e+1)[0] != delim)) ++e;
        string value = unescape_delim(extract(edit_csa, b, e));
        Edit edit;
        //cerr << b << " " << e << endl;
        edit.ParseFromString(value);
        //cerr << pb2json(edit) << endl;
        edits.push_back(edit);
    }
    return edits;
}

ostream& Counter::as_table(ostream& out, bool show_edits) {
    // write the coverage as a vector
    for (size_t i = 0; i < coverage_civ.size(); ++i) {
        out << i << "\t" << coverage_civ[i];
        if (show_edits) {
            out << "\t" << count(edit_csa, pos_key(i));
            for (auto& edit : edits_at_position(i)) out << " " << pb2json(edit);
        }
        out << endl;
    }
}

ostream& Counter::show_structure(ostream& out) {
    out << coverage_civ << endl; // graph coverage (compacted coverage_dynamic)
    out << edit_csa << endl;
    //out << " i SA ISA PSI LF BWT    T[SA[i]..SA[i]-1]" << endl;
    //csXprintf(cout, "%2I %2S %3s %3P %2p %3B   %:3T", edit_csa);
    return out;
}

}

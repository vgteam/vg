#include <thread>
#include <vg/io/protobuf_iterator.hpp>
#include "packer.hpp"
#include "statistics.hpp"
#include "../vg.hpp"

//#define debug
using namespace vg::io;

namespace vg {

const int Packer::maximum_quality = 60;
const int Packer::lru_cache_size = 4096;

size_t Packer::estimate_data_width(size_t expected_coverage) {
    return std::ceil(std::log2(2 * expected_coverage));
}

size_t Packer::estimate_batch_size(size_t num_threads) {
    size_t batch_size = max((size_t)128, (size_t)(pow(2, 14 - log2(num_threads))));
    if (batch_size % 2 != 0) {
        ++batch_size;
    }
    return batch_size;
}

size_t Packer::estimate_bin_count(size_t num_threads) {
    return pow(2, log2(num_threads) + 14);
}

Packer::Packer(const HandleGraph* graph) : graph(graph), data_width(8), cov_bin_size(0), edge_cov_bin_size(0), num_bases_dynamic(0), base_locks(nullptr), num_edges_dynamic(0), edge_locks(nullptr), node_quality_locks(nullptr), tmpfstream_locks(nullptr) { }

Packer::Packer(const HandleGraph* graph, bool record_bases, bool record_edges, bool record_edits, bool record_qualities, size_t bin_size, size_t coverage_bins, size_t data_width) :
    graph(graph), data_width(data_width), bin_size(bin_size), record_bases(record_bases), record_edges(record_edges), record_edits(record_edits), record_qualities(record_qualities) {
    // get the size of the base coverage counter
    num_bases_dynamic = 0;
    if (record_bases) {
        graph->for_each_handle([&](const handle_t& handle) { num_bases_dynamic += graph->get_length(handle); });
    }
    // get the size of the edge coverage counter
    num_edges_dynamic = 0;
    if (record_edges) {
        const VectorizableHandleGraph* vec_graph = dynamic_cast<const VectorizableHandleGraph*>(graph);
        assert(vec_graph != nullptr);      
        graph->for_each_edge([&](const edge_t& edge) {
                auto edge_index = vec_graph->edge_index(edge);
#ifdef debug
                cerr << "Observed edge at index " << edge_index << endl;
#endif
                num_edges_dynamic = std::max(num_edges_dynamic, edge_index);
            });
        ++num_edges_dynamic; // add one so our size is greater than the max element
    }
    // get the size of the node qualitiy counter
    num_nodes_dynamic = 0;
    size_t qual_coverage_bins = 0;
    if (record_qualities) {
        num_nodes_dynamic = graph->get_node_count() + 1; // add one for 1-based ranks
        if (num_bases_dynamic > 0) {
            // scale down the quality bins to be proportional in size with the base coverage bins if we can
            qual_coverage_bins = coverage_bins * ((double)num_nodes_dynamic / (double)num_bases_dynamic);
            qual_coverage_bins = max((size_t)1, qual_coverage_bins);
        }
    }
    
    // only bin if we need to
    if (num_edges_dynamic <= coverage_bins || num_bases_dynamic <= coverage_bins) {
        coverage_bins = 1;
        qual_coverage_bins = 1;
    }
    assert(coverage_bins > 0);

    // coverage counter for each bin (totally independent from the edit coverage bins)
    // they are initialized on-demand to better support sparse use-cases
    coverage_dynamic.resize(coverage_bins, nullptr);
    edge_coverage_dynamic.resize(coverage_bins, nullptr);
    node_quality_dynamic.resize(qual_coverage_bins, nullptr);
    
    // need this for every lookup, so store here
    cov_bin_size = coverage_dynamic.size() > 0 ? num_bases_dynamic / coverage_dynamic.size() : 0;
    edge_cov_bin_size = edge_coverage_dynamic.size() > 0 ? num_edges_dynamic / edge_coverage_dynamic.size() : 0;
    node_qual_bin_size = node_quality_dynamic.size() > 0 ? num_nodes_dynamic / node_quality_dynamic.size() : 0;

    // mutexes for coverage
    base_locks = new std::mutex[coverage_dynamic.size()];
    edge_locks = new std::mutex[edge_coverage_dynamic.size()];
    node_quality_locks = new std::mutex[node_quality_dynamic.size()];
    tmpfstream_locks = nullptr;
    
    // count the bins if binning
    if (bin_size) {
        n_bins = num_bases_dynamic / bin_size + 1;
    }
    if (record_edits) { 
        tmpfstream_locks = new std::mutex[n_bins];
        // open tmpfile if needed
        ensure_edit_tmpfiles_open();
    }

    // speed up quality computation if necessary
    for (size_t i = 0; i < get_thread_count(); ++i) {
        quality_cache.push_back(new LRUCache<pair<int, int>, int>(lru_cache_size));
    }
    
#ifdef debug
    cerr << "Packing across " << num_edges_dynamic << " edge slots and " << num_bases_dynamic << " base slots in " << coverage_bins << " bins" << endl;
#endif
}

void Packer::clear() {
    for (auto& counter : coverage_dynamic) {
        delete counter;
        counter = nullptr;
    }
    for (auto& counter : edge_coverage_dynamic) {
        delete counter;
        counter = nullptr;
    }
    for (auto& counter : node_quality_dynamic) {
        delete counter;
        counter = nullptr;
    }
    delete [] base_locks;
    base_locks = nullptr;
    delete [] edge_locks;
    edge_locks = nullptr;
    delete [] node_quality_locks;
    node_quality_locks = nullptr;
    delete [] tmpfstream_locks;
    tmpfstream_locks = nullptr;
    close_edit_tmpfiles();
    remove_edit_tmpfiles();
    for (auto& lru_cache : quality_cache) {
        delete lru_cache;
        lru_cache = nullptr;
    }
}

Packer::~Packer() {
    clear();
}

void Packer::load_from_file(const string& file_name) {
    ifstream in(file_name);
    if (!in) {
        stringstream ss;
        ss << "Error [Packer]: unable to read pack file: \"" << file_name << "\"" << endl;
        throw runtime_error(ss.str());
    }
    load(in);
}

void Packer::save_to_file(const string& file_name) {
    ofstream out(file_name);
    serialize(out);
}

void Packer::load(istream& in) {
    sdsl::read_member(bin_size, in);
    sdsl::read_member(n_bins, in);
    coverage_civ.load(in);
    edge_coverage_civ.load(in);
    edit_csas.resize(n_bins);
    for (size_t i = 0; i < n_bins; ++i) {
        edit_csas[i].load(in);
    }
    node_quality_civ.load(in);
    // We can only load compacted.
    is_compacted = true;
}

void Packer::merge_from_files(const vector<string>& file_names) {
#ifdef debug
    cerr << "Merging " << file_names.size() << " pack files" << endl;
#endif
    
    // load into our dynamic structures, then compact
    bool first = true;
    for (auto& file_name : file_names) {
        Packer c;
        ifstream f(file_name);
        c.load(f);
        // take bin size and counts from the first, assume they are all the same
        if (first) {
            bin_size = c.get_bin_size();
            n_bins = c.get_n_bins();
            ensure_edit_tmpfiles_open();
            first = false;
        } else {
            assert(bin_size == c.get_bin_size());
            assert(n_bins == c.get_n_bins());
        }
        c.write_edits(tmpfstreams);
        collect_coverage({&c});
    }
}

void Packer::merge_from_dynamic(vector<Packer*>& packers) {
    // load dynamic packs into our dynamic structures, then compact
    bool first = true;
    for (auto& p : packers) {
        auto& c = *p;
        c.close_edit_tmpfiles(); // flush and close temporaries
        // take bin size and counts from the first, assume they are all the same
        if (first) {
            bin_size = c.get_bin_size();
            n_bins = c.get_n_bins();
            ensure_edit_tmpfiles_open();
            first = false;
        } else {
            assert(bin_size == c.get_bin_size());
            assert(n_bins == c.get_n_bins());
        }
        c.write_edits(tmpfstreams);
    }
    collect_coverage(packers);
}

size_t Packer::get_bin_size(void) const {
    return bin_size;
}

size_t Packer::get_n_bins(void) const {
    return n_bins;
}

size_t Packer::bin_for_position(size_t i) const {
    if (bin_size > 0) {
        return i / bin_size;
    } else {
        return 0;
    }
}

void Packer::write_edits(vector<ofstream*>& out) const {
    for (size_t i = 0; i < n_bins; ++i) {
        write_edits(*out[i], i);
    }
}

void Packer::write_edits(ostream& out, size_t bin) const {
    if (is_compacted) {
        out << extract(edit_csas[bin], 0, edit_csas[bin].size()-2) << delim1; // chomp trailing null, add back delim        
    } else {
        // uncompacted, so just cat the edit file for this bin onto out
        if (edit_tmpfile_names.size()) {
            ifstream edits(edit_tmpfile_names[bin], std::ios_base::binary);
            out << edits.rdbuf() << delim1;
        }
    }
}

void Packer::collect_coverage(const vector<Packer*>& packers) {
    // assume the same basis vector
    assert(!is_compacted);
    if (record_bases) {
#pragma omp parallel for
        for (size_t i = 0; i < coverage_dynamic.size(); ++i) {
            size_t base_offset = i * cov_bin_size;
            size_t bin_size = coverage_bin_size(i);
            for (size_t j = 0; j < bin_size; ++j) {
                size_t inc_cov = 0;
                for (size_t k = 0; k < packers.size(); ++k) {
                    inc_cov += packers[k]->coverage_at_position(j + base_offset);
                }
                increment_coverage(j + base_offset, inc_cov);
            }
        }
    }
    if (record_edges) {
#pragma omp parallel for
        for (size_t i = 0; i < edge_coverage_dynamic.size(); ++i) {
            size_t edge_base_offset = i * edge_cov_bin_size;
            size_t edge_bin_size = edge_coverage_bin_size(i);
            for (size_t j = 0; j < edge_bin_size; ++j) {
                size_t inc_edge_cov = 0;
                for (size_t k = 0; k < packers.size(); ++k) {
                    inc_edge_cov += packers[k]->edge_coverage(j + edge_base_offset);
                }
                increment_edge_coverage(j + edge_base_offset, inc_edge_cov);
            }
        }
    }
    if (record_qualities) {
#pragma omp parallel for
        for (size_t i = 1; i < node_quality_dynamic.size(); ++i) {
            size_t qual_base_offset = i * node_qual_bin_size;
            size_t qual_bin_size = node_quality_bin_size(i);
            for (size_t j = 0; j < qual_bin_size; ++j) {
                size_t inc_qual_cov = 0;
                for (size_t k = 0; k < packers.size(); ++k) {
                    inc_qual_cov += packers[k]->total_node_quality(j + qual_base_offset);
                }
                increment_node_quality(j + qual_base_offset, inc_qual_cov);
            }
        }            
    }
}

size_t Packer::serialize(std::ostream& out,
                          sdsl::structure_tree_node* s,
                          std::string name) {
    make_compact();
    sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(s, name, sdsl::util::class_name(*this));
    size_t written = 0;
    written += sdsl::write_member(bin_size, out, child, "bin_size_" + name);
    written += sdsl::write_member(edit_csas.size(), out, child, "n_bins_" + name);
    written += coverage_civ.serialize(out, child, "graph_coverage_" + name);
    written += edge_coverage_civ.serialize(out, child, "edge_coverage_" +name);
    for (auto& edit_csa : edit_csas) {
        written += edit_csa.serialize(out, child, "edit_csa_" + name);
    }
    written += node_quality_civ.serialize(out, child, "node_quality_" + name);
    sdsl::structure_tree::add_size(child, written);
    return written;
}

void Packer::make_compact(void) {
    // pack the dynamic countarray and edit coverage into the compact data structure
    if (is_compacted) {
#ifdef debug
        cerr << "Packer is already compact" << endl;
#endif
        return;
    } else {
#ifdef debug
        cerr << "Need to make packer compact" << endl;
#endif
    }
    // sync edit file
    close_edit_tmpfiles();
    
    // temporaries for construction
    size_t basis_length = coverage_size();
    int_vector<> coverage_iv;
    size_t edge_coverage_length = edge_vector_size();
    int_vector<> edge_coverage_iv;
    size_t node_quality_length = node_quality_vector_size();
    int_vector<> node_quality_iv;
    
#ifdef debug
    cerr << "Concatenating entries for " << basis_length << " bases and " << edge_coverage_length << " edges" << endl;
#endif
    
#pragma omp parallel
    {
#pragma omp single
        {
#pragma omp task
            {
                util::assign(coverage_iv, int_vector<>(basis_length));
            }
#pragma omp task
            {
                util::assign(edge_coverage_iv, int_vector<>(edge_coverage_length));
            }
#pragma omp task
            {
                util::assign(node_quality_iv, int_vector<>(node_quality_length));
            }
        }
    }
#pragma omp parallel for
    for (size_t i = 0; i < basis_length; ++i) {
        coverage_iv[i] = coverage_at_position(i);
    }
#pragma omp parallel for
    for (size_t i = 0; i < edge_coverage_length; ++i) {
        edge_coverage_iv[i] = edge_coverage(i);
    }
#pragma omp parallel for
    for (size_t i = 1; i < node_quality_length; ++i) {
        node_quality_iv[i] = average_node_quality(i);
    }

    #pragma omp parallel
    {
#pragma omp single
        {
#pragma omp task
            {
                util::assign(coverage_civ, coverage_iv);
            }
#pragma omp task
            {
                util::assign(edge_coverage_civ, edge_coverage_iv);
            }
#pragma omp task
            {
                util::assign(node_quality_civ, node_quality_iv);
            }
        }
    }
    
    edit_csas.resize(edit_tmpfile_names.size());
    construct_config::byte_algo_sa = SE_SAIS;
#pragma omp parallel for
    for (size_t i = 0; i < edit_tmpfile_names.size(); ++i) {
        construct(edit_csas[i], edit_tmpfile_names[i], 1);
    }
    // construct the record marker bitvector
    remove_edit_tmpfiles();
    is_compacted = true;
}

void Packer::make_dynamic(void) {
    if (!is_compacted) return;
    // unpack the compact represenation into the countarray
    assert(false); // not implemented
    is_compacted = false;
}

bool Packer::is_dynamic(void) const {
    return !is_compacted;
}

const HandleGraph* Packer::get_graph() const {
    return graph;
}

void Packer::ensure_edit_tmpfiles_open(void) {
    if (tmpfstreams.empty()) {
        string base = "vg-pack_";
        string edit_tmpfile_name = temp_file::create(base);
        temp_file::remove(edit_tmpfile_name); // remove this; we'll use it as a base name
        // for as many bins as we have, make a temp file
        tmpfstreams.resize(n_bins);
        edit_tmpfile_names.resize(n_bins);
        for (size_t i = 0; i < n_bins; ++i) {
            edit_tmpfile_names[i] = edit_tmpfile_name+"_"+convert(i);
            tmpfstreams[i] = new ofstream;
            tmpfstreams[i]->open(edit_tmpfile_names[i], std::ios_base::binary);
            assert(tmpfstreams[i]->is_open());
        }
    }
}

void Packer::close_edit_tmpfiles(void) {
    if (!tmpfstreams.empty()) {
        for (auto& tmpfstream : tmpfstreams) {
            *tmpfstream << delim1; // pad
            tmpfstream->close();
            delete tmpfstream;
        }
        tmpfstreams.clear();
    }
}

void Packer::remove_edit_tmpfiles(void) {
    if (!edit_tmpfile_names.empty()) {
        for (auto& name : edit_tmpfile_names) {
            std::remove(name.c_str());
        }
        edit_tmpfile_names.clear();
    }
}

void Packer::add(const Alignment& aln, int min_mapq, int min_baseq, int trim_ends) {
    // mapping quality threshold filter
    int mapping_quality = aln.mapping_quality();
    if (mapping_quality < min_mapq) {
        return;
    }
    // count the nodes, edges, and edits
    Mapping prev_mapping;
    bool has_prev_mapping = false;
    int prev_bq_total = 0;
    int prev_bq_count = 0;
    size_t position_in_read = 0;

    size_t trim_start_mapping = 0;
    size_t trim_start_edit = 0;
    size_t trim_start_offset = 0;
    size_t trim_end_mapping = numeric_limits<size_t>::max();
    size_t trim_end_edit = numeric_limits<size_t>::max();
    size_t trim_end_offset = numeric_limits<size_t>::max();
    if (trim_ends > 0) {
        // do a pass to determine trimming.
        // we count "trim_ends" aligned bases from each end, and remember
        // the mapping/edit/offset when the count is reached.
        // in the logic below, everything outside these intervals will
        // be ignored.  
        size_t ab = 0;
        size_t found = false;
        // todo: we could do this pass in the main loop down below...
        for (size_t i = 0; i < aln.path().mapping_size() && !found; ++i) {
            const Mapping& mapping = aln.path().mapping(i);
            for (size_t j = 0; j < mapping.edit_size() && !found; ++j) {
                const Edit& edit = mapping.edit(j);
                if (edit_is_match(edit)) {
                    if (ab + edit.to_length() >= trim_ends) {
                        trim_start_mapping = i;
                        trim_start_edit = j;
                        trim_start_offset = trim_ends - ab;
                        found = true;
                    }
                    ab += edit.to_length();
                }
            }            
        }
        // same as above, but backwards
        ab = 0;
        found = false;
        for (size_t i = 0; i < aln.path().mapping_size() && !found; ++i) {
            const Mapping& mapping = aln.path().mapping(aln.path().mapping_size() - 1 - i);
            for (size_t j = 0; j < mapping.edit_size() && !found; ++j) {
                const Edit& edit = mapping.edit(mapping.edit_size() - 1 - j);
                if (edit_is_match(edit)) {
                    if (ab + edit.to_length() >= trim_ends) {
                        trim_end_mapping = aln.path().mapping_size() - 1 - i;
                        trim_end_edit = mapping.edit_size() - 1 - j;
                        trim_end_offset = edit.to_length() - (trim_ends - ab);
                        found = true;
                    }
                    ab += edit.to_length();
                }
            }            
        }
    }
    function<bool(size_t, size_t, size_t)> trim_check = [&](size_t mi, size_t ei, size_t eo) {
        if (mi > trim_start_mapping ||
            (mi == trim_start_mapping && ei > trim_start_edit) ||
            (mi == trim_start_mapping && ei == trim_start_edit && eo >= trim_start_offset)) {
            if (mi < trim_end_mapping ||
                (mi == trim_end_mapping && ei < trim_end_edit) ||
                (mi == trim_end_mapping && ei == trim_end_edit && eo < trim_end_offset)) {
                return true;
            }
        }
        return false;
    };

    size_t cur_pos = 0;
    for (size_t mi = 0; mi < aln.path().mapping_size(); ++mi) {
        auto& mapping = aln.path().mapping(mi);
        if (!mapping.has_position()) {
#ifdef debug
            cerr << "Mapping has no position" << endl;
#endif
            has_prev_mapping = false;
            continue;
        }
        // skip nodes outside of our graph, assuming this may be a subgraph
        if (!graph->has_node(mapping.position().node_id())) {
            has_prev_mapping = false;
            continue;
        }
        size_t i = position_in_basis(mapping.position());
        size_t node_quality_index = node_index(mapping.position().node_id());
        size_t total_node_quality = 0;
        // keep track of average base quality in the mapping
        int bq_total = 0;
        int bq_count = 0;
        int ei = 0;
        if (record_bases || record_qualities) {
            for (auto& edit : mapping.edit()) {
                if (edit_is_match(edit)) {
#ifdef debug
                    cerr << "Recording a match" << endl;
#endif
                    int direction = mapping.position().is_reverse() ? -1 : 1;
                    for (size_t j = 0; j < edit.from_length(); ++j, ++position_in_read) {
                        int64_t coverage_idx = i + direction * j;
                        int base_quality = compute_quality(aln, position_in_read);
                        bq_total += base_quality;
                        ++bq_count;
                        // base quality threshold filter (only if we found some kind of quality)
                        if (record_bases && (base_quality < 0 || base_quality >= min_baseq) && (trim_ends <= 0 || trim_check(mi, ei, j))) {
                            increment_coverage(coverage_idx);
                            if (record_qualities && mapping_quality > 0) {
                                total_node_quality += mapping_quality;
                            }
                        }
                    }         
                } else if (record_edits) {
                    // we represent things on the forward strand
                    string pos_repr = pos_key(i);
                    string edit_repr = edit_value(edit, mapping.position().is_reverse());
                    size_t bin = bin_for_position(i);
                    std::lock_guard<std::mutex> guard(tmpfstream_locks[bin]);
                    *tmpfstreams[bin] << pos_repr << edit_repr;
                } 
                if (mapping.position().is_reverse()) {
                    i -= edit.from_length();
                } else {
                    i += edit.from_length();
                }
                if (!edit_is_match(edit)) {
                    position_in_read += edit.to_length();
                }
                ++ei;
            }
            if (total_node_quality > 0) {
                increment_node_quality(node_quality_index, total_node_quality);
            }
        }
        
        if (record_edges && has_prev_mapping && prev_mapping.position().node_id() != mapping.position().node_id() &&
            (trim_ends <= 0 || (trim_check(mi, 0, 0) &&
                                trim_check(mi - 1, aln.path().mapping(mi - 1).edit_size() - 1, aln.path().mapping(mi - 1).edit(aln.path().mapping(mi - 1).edit_size() - 1).to_length() - 1)))) {
            // Note: we are effectively ignoring edits here.  So an edge is covered even
            // if there's a sub or indel at either of its ends in the path.  
            Edge e;
            e.set_from(prev_mapping.position().node_id());
            e.set_from_start(prev_mapping.position().is_reverse());
            e.set_to(mapping.position().node_id());
            e.set_to_end(mapping.position().is_reverse());
            size_t edge_idx = edge_index(e);
#ifdef debug
            cerr << "Observed visit to edge " << pb2json(e) << " at index " << edge_idx << endl;
#endif
            if (edge_idx != 0) {
                // heuristic:  for an edge, we average out the base qualities from the matches in its two flanking mappings
                int avg_base_quality = -1;
                if (!aln.quality().empty()) {
                    if (bq_count + prev_bq_count == 0) {
                        avg_base_quality = 0;
                    } else {
                        avg_base_quality = (float)(bq_total + prev_bq_total) / (bq_count + prev_bq_count);
                    }
                }
                // base quality threshold filter (only if we found some kind of quality)
                if (avg_base_quality < 0 || avg_base_quality >= min_baseq) {
                    increment_edge_coverage(edge_idx);
                }
            }
        }

        prev_mapping = mapping;
        has_prev_mapping = true;
    }
}

// find the position on the forward strand in the sequence vector
size_t Packer::position_in_basis(const Position& pos) const {
    // get position on the forward strand
    if (pos.is_reverse()) {
        return (int64_t)dynamic_cast<const VectorizableHandleGraph*>(graph)->node_vector_offset(pos.node_id())
            + (int64_t)reverse(pos, graph->get_length(graph->get_handle(pos.node_id()))).offset() - 1;
    } else {
        return (int64_t)dynamic_cast<const VectorizableHandleGraph*>(graph)->node_vector_offset(pos.node_id())
            + (int64_t)pos.offset();
    }
}

string Packer::pos_key(size_t i) const {
    Position pos;
    size_t offset = 2;
    pos.set_node_id(i+offset);
    string pos_repr;
    pos.SerializeToString(&pos_repr);
    stringstream s;
    s << delim1 << delim2 << delim1 << escape_delims(pos_repr);
    return s.str();
}

string Packer::edit_value(const Edit& edit, bool revcomp) const {
    string edit_repr;
    if (revcomp) {
        reverse_complement_edit(edit).SerializeToString(&edit_repr);
    } else {
        edit.SerializeToString(&edit_repr);
    }
    stringstream s;
    s << delim1 << escape_delims(edit_repr);
    return s.str();
}

string Packer::escape_delims(const string& s) const {
    return escape_delim(escape_delim(s, delim1), delim2);
}

string Packer::unescape_delims(const string& s) const {
    return unescape_delim(unescape_delim(s, delim1), delim2);
}

string Packer::escape_delim(const string& s, char d) const {
    string escaped; escaped.reserve(s.size());
    for (size_t i = 0; i < s.size(); ++i) {
        char c = s[i];
        escaped.push_back(c);
        if (c == d) escaped.push_back(c);
    }
    return escaped;
}

string Packer::unescape_delim(const string& s, char d) const {
    string unescaped; unescaped.reserve(s.size());
    for (size_t i = 0; i < s.size()-1; ++i) {
        char c = s[i];
        char b = s[i+1];
        if (c == d && b == d) {
            unescaped.push_back(c);
        } else {
            unescaped.push_back(c);
            if (i == s.size()-2) unescaped.push_back(b);
        }
    }
    return unescaped;
}

size_t Packer::coverage_size(void) const {
    if (is_compacted){
        return coverage_civ.size();
    }
    else{
        return num_bases_dynamic;
    }
}

size_t Packer::edge_vector_size(void) const{
    if (is_compacted){
        return edge_coverage_civ.size();
    }
    else{
        return num_edges_dynamic;
    }
}

size_t Packer::node_quality_vector_size(void) const {
    if (is_compacted) {
        return node_quality_civ.size();
    } else {
        return num_nodes_dynamic;
    }
}

pair<size_t, size_t> Packer::coverage_bin_offset(size_t i) const {
    size_t bin = min((size_t)(i / cov_bin_size), (size_t)(coverage_dynamic.size() - 1));
    // last bin can have different size so we don't use mod
    size_t offset = i - bin * cov_bin_size;
    return make_pair(bin, offset);
}

pair<size_t, size_t> Packer::edge_coverage_bin_offset(size_t i) const {
    size_t bin = min((size_t)(i / edge_cov_bin_size), (size_t)(edge_coverage_dynamic.size() - 1));
    // last bin can have different size so we don't use mod
    size_t offset = i - bin * edge_cov_bin_size;
    return make_pair(bin, offset);
}

pair<size_t, size_t> Packer::node_quality_bin_offset(size_t i) const {
    size_t bin = min((size_t)(i / node_qual_bin_size), (size_t)(node_quality_dynamic.size() - 1));
    // last bin can have different size so we don't use mod
    size_t offset = i - bin * node_qual_bin_size;
    return make_pair(bin, offset);
}

size_t Packer::coverage_bin_size(size_t i) const {
    size_t bin_size = cov_bin_size;
    if (i == coverage_dynamic.size() - 1) {
        bin_size += num_bases_dynamic % coverage_dynamic.size();
    }
    return bin_size;
}

size_t Packer::edge_coverage_bin_size(size_t i) const {
    size_t bin_size = edge_cov_bin_size;
    if (i == edge_coverage_dynamic.size() - 1) {
        bin_size += num_edges_dynamic % edge_coverage_dynamic.size();
    }
    return bin_size;
}

size_t Packer::node_quality_bin_size(size_t i) const {
    size_t bin_size = node_qual_bin_size;
    if (i == node_quality_dynamic.size() - 1) {
        bin_size += num_nodes_dynamic % node_quality_dynamic.size();
    }
    return bin_size;
}

void Packer::init_coverage_bin(size_t i) {
    if (coverage_dynamic[i] == nullptr) {
        coverage_dynamic[i] = new gcsa::CounterArray(coverage_bin_size(i), data_width);
    }
}

void Packer::init_edge_coverage_bin(size_t i) {
    if (edge_coverage_dynamic[i] == nullptr) {
        edge_coverage_dynamic[i] = new gcsa::CounterArray(edge_coverage_bin_size(i), data_width);
    }
}

void Packer::init_node_quality_bin(size_t i) {
    if (node_quality_dynamic[i] == nullptr) {
        // add some bits to the data width because we have to add a quality and not 1 each time
        node_quality_dynamic[i] = new gcsa::CounterArray(node_quality_bin_size(i), data_width + 12);
    }
}

void Packer::increment_coverage(size_t i) {
    pair<size_t, size_t> bin_offset = coverage_bin_offset(i);
    std::lock_guard<std::mutex> guard(base_locks[bin_offset.first]);
    init_coverage_bin(bin_offset.first);
    coverage_dynamic.at(bin_offset.first)->increment(bin_offset.second);
}

void Packer::increment_coverage(size_t i, size_t v) {
    if (v > 0) {
        pair<size_t, size_t> bin_offset = coverage_bin_offset(i);
        std::lock_guard<std::mutex> guard(base_locks[bin_offset.first]);
        init_coverage_bin(bin_offset.first);
        coverage_dynamic.at(bin_offset.first)->increment(bin_offset.second, v);
    }
}

void Packer::increment_edge_coverage(size_t i) {
    pair<size_t, size_t> bin_offset = edge_coverage_bin_offset(i);
    std::lock_guard<std::mutex> guard(edge_locks[bin_offset.first]);
    init_edge_coverage_bin(bin_offset.first);
    edge_coverage_dynamic.at(bin_offset.first)->increment(bin_offset.second);
#ifdef debug
    cerr << "Set coverage of edge " << i << " at " << bin_offset.first << "/" << bin_offset.second << " to " << (*edge_coverage_dynamic.at(bin_offset.first))[bin_offset.second] << endl;
#endif
}

void Packer::increment_edge_coverage(size_t i, size_t v) {
    if (v > 0) {
        pair<size_t, size_t> bin_offset = edge_coverage_bin_offset(i);
        std::lock_guard<std::mutex> guard(edge_locks[bin_offset.first]);
        init_edge_coverage_bin(bin_offset.first);
        edge_coverage_dynamic.at(bin_offset.first)->increment(bin_offset.second, v);
#ifdef debug
        cerr << "Set coverage of edge " << i << " at " << bin_offset.first << "/" << bin_offset.second << " to " << (*edge_coverage_dynamic.at(bin_offset.first))[bin_offset.second] << endl;
#endif
    }
}

void Packer::increment_node_quality(size_t i, size_t v) {
    if (v > 0) {
        pair<size_t, size_t> bin_offset = node_quality_bin_offset(i);
        std::lock_guard<std::mutex> guard(node_quality_locks[bin_offset.first]);
        init_node_quality_bin(bin_offset.first);
        node_quality_dynamic.at(bin_offset.first)->increment(bin_offset.second, v);
#ifdef debug
        cerr << "Set quality of node " << i << " at " << bin_offset.first << "/" << bin_offset.second << " to " << (*node_quality_dynamic.at(bin_offset.first))[bin_offset.second] << endl;
#endif
    }        
}

bool Packer::has_qualities() const {
    if (is_compacted) {
        for (size_t i = 0; i < node_quality_civ.size(); ++i) {
            if (node_quality_civ[i] > 0) {
                return true;
            }
        }
    } else {
        for (size_t i = 0; i < node_quality_dynamic.size(); ++i) {
            if (node_quality_dynamic[i] != nullptr) {
                for (size_t j = 0; j < node_quality_dynamic[i]->size(); ++j) {
                    if ((*node_quality_dynamic.at(i))[j] > 0) {
                        return true;
                    }
                }
            }
        }
    }
    return false;
}

size_t Packer::coverage_at_position(size_t i) const {
    if (is_compacted) {
        return coverage_civ[i];
    } else {
        pair<size_t, size_t> bin_offset = coverage_bin_offset(i);
        if (coverage_dynamic[bin_offset.first] == nullptr) {
            return 0;
        } else {
            return (*coverage_dynamic.at(bin_offset.first))[bin_offset.second];
        }
    }
}

size_t Packer::edge_coverage(size_t i) const {
    if (is_compacted){
        return edge_coverage_civ[i];
    }
    else{
        pair<size_t, size_t> bin_offset = edge_coverage_bin_offset(i);
        if (edge_coverage_dynamic[bin_offset.first] == nullptr) {
            return 0;
        } else {
            return (*edge_coverage_dynamic.at(bin_offset.first))[bin_offset.second];
        }
    }
}
            
size_t Packer::edge_coverage(Edge& e) const {
    size_t pos = edge_index(e);
    return edge_coverage(pos);
}

size_t Packer::total_node_quality(size_t i) const {
    if (is_compacted) {
        size_t avg_qual = average_node_quality(i);
        Position pos;
        pos.set_node_id(index_to_node(i));
        pos.set_is_reverse(false);
        size_t cov_pos = position_in_basis(pos);
        size_t node_len = graph->get_length(graph->get_handle(pos.node_id()));
        size_t base = position_in_basis(pos);
        size_t coverage = 0;
        for (size_t i = 0; i < node_len; ++i) {
            coverage += coverage_at_position(base + i);
        }
        return avg_qual * coverage;
    } else {
        pair<size_t, size_t> bin_offset = node_quality_bin_offset(i);
        if (node_quality_dynamic[bin_offset.first] == nullptr) {
            return 0;
        } else {
            return (*node_quality_dynamic.at(bin_offset.first))[bin_offset.second];
        }
    }
}

size_t Packer::average_node_quality(size_t i) const {
    if (is_compacted) {
        return node_quality_civ[i];
    } else {
        Position pos;
        pos.set_node_id(index_to_node(i));
        pos.set_is_reverse(false);
        size_t cov_pos = position_in_basis(pos);
        size_t node_len = graph->get_length(graph->get_handle(pos.node_id()));
        size_t base = position_in_basis(pos);
        size_t total_coverage = 0;
        for (size_t i = 0; i < node_len; ++i) {
            total_coverage += coverage_at_position(base + i);
        }
        if (total_coverage == 0) {
            assert(total_node_quality(i) == 0);
            return 0;
        } else {
            return (size_t)std::lround(total_node_quality(i) / total_coverage);
        }
    }
}

vector<Edit> Packer::edits_at_position(size_t i) const {
    vector<Edit> edits;
    if (i == 0) return edits;
    string key = pos_key(i);
    size_t bin = bin_for_position(i);
    auto& edit_csa = edit_csas[bin];
    auto occs = locate(edit_csa, key);
    for (size_t i = 0; i < occs.size(); ++i) {
        // walk from after the key and delim1 to the next end-sep
        size_t b = occs[i] + key.size() + 1;
        size_t e = b;
        // look for an odd number of delims
        // run until we find a delim
        while (true) {
            while (extract(edit_csa, e, e)[0] != delim1) ++e;
            // now we are matching the delim... count them
            size_t f = e;
            while (extract(edit_csa, f, f)[0] == delim1) ++f;
            size_t c = f - e;
            e = f; // set pointer to last delim
            if (c % 2 != 0) {
                break;
            }
        }
        string value = unescape_delims(extract(edit_csa, b, e));
        Edit edit;
        vg::io::ProtobufIterator<Edit>::parse_from_string(edit, value);
        edits.push_back(edit);
    }
    return edits;
}

size_t Packer::edge_index(const Edge& e) const {
    edge_t edge = graph->edge_handle(graph->get_handle(e.from(), e.from_start()),
                                     graph->get_handle(e.to(), e.to_end()));
                                     
    if (!graph->has_edge(edge)) {
        // We can only query the edge index for edges that exist. This edge doesn't.
        return 0;
    }
                                     
    return dynamic_cast<const VectorizableHandleGraph*>(graph)->edge_index(edge);
}

size_t Packer::node_index(nid_t node_id) const {
    if (!graph->has_node(node_id)) {
        return 0;
    }
    return dynamic_cast<const VectorizableHandleGraph*>(graph)->id_to_rank(node_id);
}

nid_t Packer::index_to_node(size_t i) const {
    return dynamic_cast<const VectorizableHandleGraph*>(graph)->rank_to_id(i);
}

ostream& Packer::as_table(ostream& out, bool show_edits, vector<vg::id_t> node_ids) {
#ifdef debug
    cerr << "Packer table of " << coverage_civ.size() << " rows:" << endl;
#endif

    out << "seq.pos" << "\t"
        << "node.id" << "\t"
        << "node.offset" << "\t"
        << "coverage";
    if (show_edits) out << "\t" << "edits";
    out << endl;
    // write the coverage as a vector
    for (size_t i = 0; i < coverage_civ.size(); ++i) {
        nid_t node_id = dynamic_cast<const VectorizableHandleGraph*>(graph)->node_at_vector_offset(i+1);
        if (!node_ids.empty() && find(node_ids.begin(), node_ids.end(), node_id) == node_ids.end()) {
            continue;
        }
        size_t offset = i - dynamic_cast<const VectorizableHandleGraph*>(graph)->node_vector_offset(node_id);
        out << i << "\t" << node_id << "\t" << offset << "\t" << coverage_civ[i];
        if (show_edits) {
            out << "\t" << count(edit_csas[bin_for_position(i)], pos_key(i));
            for (auto& edit : edits_at_position(i)) out << " " << pb2json(edit);
        }
        out << endl;
    }
    return out;
}

ostream& Packer::as_edge_table(ostream& out, vector<vg::id_t> node_ids) {
#ifdef debug
    cerr << "Packer edge table of " << edge_coverage_civ.size() << " rows:" << endl;
#endif

    out << "from.id" << "\t"
        << "from.start" << "\t"
        << "to.id" << "\t"
        << "to.end" << "\t"
        << "coverage" << endl;
    graph->for_each_edge([&](const edge_t& handle_edge) {
            Edge edge;
            edge.set_from(graph->get_id(handle_edge.first));
            edge.set_from_start(graph->get_is_reverse(handle_edge.first));
            edge.set_to(graph->get_id(handle_edge.second));
            edge.set_to_end(graph->get_is_reverse(handle_edge.second));
            
            if (!node_ids.empty() &&
                (find(node_ids.begin(), node_ids.end(), edge.from()) == node_ids.end() ||
                find(node_ids.begin(), node_ids.end(), edge.to()) == node_ids.end())) {
                
                // We need to skip this edge because it deals with nodes outsode of our set.
                return true;
            }
            
            // Otherwise, we need to use the edge; all edges are visited exactly once.
            // But we want to output it smaller node ID first, which is not guaranteed.
            if (edge.from() > edge.to()) {
                {
                    // Swap the from and to around
                    nid_t temp = edge.from();
                    edge.set_from(edge.to());
                    edge.set_to(temp);
                }
                {
                    // And the flags 
                    bool temp = edge.from_start();
                    edge.set_from_start(!edge.to_end());
                    edge.set_to_end(!temp);
                }
            }
            
            // TODO: we don't canonicalize self loops at all
            
            // Print out the edge
            out << edge.from() << "\t"
                << edge.from_start() << "\t"
                << edge.to() << "\t"
                << edge.to_end() << "\t"
                << edge_coverage_civ[edge_index(edge)]
                << endl;
            
            // Look at the enxt edge
            return true;
        });
    return out;
}
    
ostream& Packer::as_quality_table(ostream& out, vector<vg::id_t> node_ids) {
#ifdef debug
    cerr << "Packer quality table of " << node_quality_civ.size() << " rows:" << endl;
#endif

    out << "node.rank" << "\t"
        << "node.id" << "\t"
        << "avg-mapq";
    out << endl;
    for (size_t i = 1; i < node_quality_civ.size(); ++i) {
        nid_t node_id = index_to_node(i);
        if (!node_ids.empty() && find(node_ids.begin(), node_ids.end(), node_id) == node_ids.end()) {
            continue;
        }
        out << i << "\t" << node_id << "\t" << node_quality_civ[i] << endl;
    }
    return out;
}

ostream& Packer::show_structure(ostream& out) {
    out << coverage_civ << endl; // graph coverage (compacted coverage_dynamic)
    for (auto& edit_csa : edit_csas) {
        out << edit_csa << endl;
    }
    //out << " i SA ISA PSI LF BWT    T[SA[i]..SA[i]-1]" << endl;
    //csXprintf(cout, "%2I %2S %3s %3P %2p %3B   %:3T", edit_csa);
    return out;
}

int Packer::compute_quality(const Alignment& aln, size_t position_in_read) const {
    int map_quality = (int)aln.mapping_quality();
    int base_quality = -1;
    if (!aln.quality().empty()) {
        base_quality = (int)aln.quality()[position_in_read];
    }
    return combine_qualities(map_quality, base_quality);
}

int Packer::combine_qualities(int map_quality, int base_quality) const {
    if (base_quality < 0) {
        // no base quality in read: just return the mapping quality
        return map_quality;
    } else {
        if (base_quality == 0 || map_quality == 0) {
            return 0;
        }

        // look up the mapping and base quality in the cache to avoid recomputing
        auto& qual_cache = *quality_cache[omp_get_thread_num()];
        pair<int, bool> cached = qual_cache.retrieve(make_pair(map_quality, base_quality));
        if (cached.second == true) {
            return cached.first;
        } else {
            // assume independence: P[Correct] = P[Correct Base] * P[Correct Map]
            // --> P[Error] = 1 - (1 - P[Base Error]) * (1 - P[Map Error])
            double p_err = logprob_invert(logprob_invert(phred_to_logprob(base_quality)) +
                                          logprob_invert(phred_to_logprob(map_quality)));
            // clamp our quality to 60
            int qual = min((int)logprob_to_phred(p_err), (int)maximum_quality);
            // update the cache
            qual_cache.put(make_pair(map_quality, base_quality), qual);
            return qual;
        }
    }
}

    
}

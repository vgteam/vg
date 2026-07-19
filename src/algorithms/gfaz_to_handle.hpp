#ifndef VG_ALGORITHMS_GFAZ_TO_HANDLE_HPP_INCLUDED
#define VG_ALGORITHMS_GFAZ_TO_HANDLE_HPP_INCLUDED

/**
 * \file gfaz_to_handle.hpp
 *
 * Defines algorithms for copying data from GFAz files into handle graphs.
 */

#include <GFAz/compress/io/gfa_decoder.hpp>

#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

#include "gfa_to_handle.hpp"

namespace vg {
namespace algorithms {
using namespace std;

/// GFAz parser implementation that adapts GFAz's decoded record visitor to
/// vg's shared GFA listener interface.
class GFAzParser : public GFAParser, private gfaz::GfaRecordVisitor {
public:
    void parse(istream& in) override;
    void parse(const string& filename) override;

    /// Return true if the named file begins with the GFAz magic number.
    /// Returns false for stdin ("-") or empty filenames.
    static bool looks_like_gfaz(const string& filename);

    /// Return true if the first sizeof(uint32_t) bytes of the buffer encode the
    /// GFAz magic number.
    static bool has_magic(const char* buffer, size_t length);

private:
    struct PendingRGFAVisit {
        nid_t id;
        int64_t offset;
        size_t length;
        string path_name;
        int64_t path_rank;
    };

    vector<PendingRGFAVisit> pending_rgfa_visits;
    unordered_map<string, int64_t> rgfa_path_ranks;

    void prepare_decode();
    void finish_decode();

    void on_header(const string& header_line) override;
    void on_segment(uint32_t id, const string& sequence,
                    const gfaz::GfaTagList& tags) override;
    void on_link(uint32_t from_id, bool from_is_reverse,
                 uint32_t to_id, bool to_is_reverse,
                 const string& overlap,
                 const gfaz::GfaTagList& tags) override;
    void on_path(const string& name, const vector<gfaz::NodeId>& visits,
                 const string& overlap,
                 const gfaz::GfaTagList& tags) override;
    void on_walk(const string& sample_name, uint32_t haplotype,
                 const string& contig_name, int64_t sequence_start,
                 int64_t sequence_end,
                 const vector<gfaz::NodeId>& visits) override;

    static visit_source_t make_visit_source(const vector<gfaz::NodeId>& visits);

    /// Wrap a function in a handler for GFADuplicatePathError that reports
    /// that a line of the given type is being skipped.
    void handle_duplicate_paths(char line_type,
                                const function<void(void)>& callback);
};

/// Load either a GFA or a GFAz into a HandleGraph.
/// TODO: Fix https://github.com/vgteam/vg/issues/4880 and remove this!
void load_gfa_or_gfaz_to_handle_graph(const string& filename,
                                      MutableHandleGraph* graph,
                                      GFAIDMapInfo* translation = nullptr);

/// Load either a GFA or a GFAz into a HandleGraph, saving the translation to the given file.
/// TODO: Fix https://github.com/vgteam/vg/issues/4880 and remove this!
void load_gfa_or_gfaz_to_handle_graph(const string& filename,
                                      MutableHandleGraph* graph,
                                      const string& translation_filename);

/// Load a GFA or GFAz into a PathHandleGraph.
/// TODO: Fix https://github.com/vgteam/vg/issues/4880 and remove this!
void load_gfa_or_gfaz_to_path_handle_graph(const string& filename,
                                           MutablePathMutableHandleGraph* graph,
                                           GFAIDMapInfo* translation = nullptr,
                                           int64_t max_rgfa_rank = numeric_limits<int64_t>::max(),
                                           unordered_set<PathSense>* ignore_sense = nullptr);

/// Load a GFA or GFAz into a PathHandleGraph, saving the translation to the given file.
/// TODO: Fix https://github.com/vgteam/vg/issues/4880 and remove this!
void load_gfa_or_gfaz_to_path_handle_graph(const string& filename,
                                           MutablePathMutableHandleGraph* graph,
                                           int64_t max_rgfa_rank,
                                           const string& translation_filename,
                                           unordered_set<PathSense>* ignore_sense = nullptr);

}
}

#endif

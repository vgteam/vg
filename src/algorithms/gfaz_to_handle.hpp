#ifndef VG_ALGORITHMS_GFAZ_TO_HANDLE_HPP_INCLUDED
#define VG_ALGORITHMS_GFAZ_TO_HANDLE_HPP_INCLUDED

/**
 * \file gfa_to_handle.hpp
 *
 * Defines algorithms for copying data from GFA files into handle graphs
 */

#include <iostream>
#include <string>

#include "gfa_to_handle.hpp"

namespace vg {
namespace algorithms {
using namespace std;

/// GFAZ parser implementation that emits the same listener events from
/// decompressed GFAZ data.
class GFAzParser : public GFAParser {
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
    /// Wrap a function in a handler for GFADuplicatePathError that reports
    /// that a line of the given line tyoe is being skipped.
    void handle_duplicate_paths(char line_type, const std::function<void(void)>& callback);
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

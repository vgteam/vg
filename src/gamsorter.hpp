#ifndef VG_GAMSORTER_HPP_INCLUDED
#define VG_GAMSORTER_HPP_INCLUDED

#include "vg.pb.h"
#include "stream.hpp"
#include "types.hpp"
#include "progressive.hpp"
#include <string>
#include <queue>
#include <sstream>
#include <functional>
#include <algorithm>
#include <iostream>
#include <set>
#include <vector>
#include <unordered_map>
#include <tuple>

/**
 * \file gamsorter.hpp
 * GAM sorting tools.
 */
using namespace std;
namespace vg {


// We need to know about the GAMIndex, but we don't actually need to hold one.
// So pre-declare it here.
class GAMIndex;

/// Provides the ability to sort a GAM, either "dumbly" (in memory), or
/// "streaming" into temporary files. Paired alignments are not necessarily
/// going to end up next to each other, so if sorting by position make sure to
/// set the position cross-references first if you want to be able to find
/// them.
class GAMSorter : public Progressive {
public:

    //////////////////
    // Main entry points
    //////////////////
    
    /// Create a GAM sorter, showing sort progress on standard error if show_progress is true.
    GAMSorter(bool show_progress = false);
    
    /// Sort a stream of GAM-format data, using temporary files, limiting the
    /// number of simultaneously open input files and the size of in-memory
    /// data. Optionally index the sorted GAM file into the given GAMIndex.
    void stream_sort(istream& gam_in, ostream& gam_out, GAMIndex* index_to = nullptr);
    
    /// Sort a stream of GAM-format data, loading it all into memory and doing
    /// a single giant sort operation.
    /// Optionally index the sorted GAM file into the given GAMIndex.
    void dumb_sort(istream& gam_in, ostream& gam_out, GAMIndex* index_to = nullptr);
    
    /// Sort a seekable input stream by doing one pass to load all the
    /// positions, sorting all the positions in memory, and doing another pass
    /// of jumping around to re-order all the reads.
    /// Optionally index the sorted GAM file into the given GAMIndex.
    void benedict_sort(istream& gam_in, ostream& gam_out, GAMIndex* index_to = nullptr);
    
    //////////////////
    // Supporting API
    //////////////////

    /// Sort a vector of alignments, in place.
    void sort(vector<Alignment>& alns) const;

    /// Return true if out of Alignments a and b, alignment a must come before alignment b, and false otherwise.
    bool less_than(const Alignment& a, const Alignment& b) const;
    
    /// Determine the minumum Position visited by an Alignment. The minumum
    /// Position is the lowest node ID visited by the alignment, with the
    /// lowest offset visited on that node ID as the offset, and the
    /// orientation set to false if the forward strand is visited, and true if
    /// only the reverse strand is visited.
    Position get_min_position(const Alignment& aln) const;

    /// Determine the minimum position visited by a Path, as for an Alignment.
    Position get_min_position(const Path& path) const;

    /// Return True if the given Position values are equal, and false otherwise.
    bool equal_to(const Position& a, const Position& b) const;

    /// Return True if position A is less than position B in our sort, and false otherwise.
    /// Position order is defined first by node ID, then by strand (forward first), and then by offset within the strand.
    /// We can't sort by actual base on the forward strand, because we need to be able to sort without knowing the graph's node lengths.
    bool less_than(const Position& a, const Position& b) const;
    
    /// Return true if out of pos_t items a and b, a must come before b, and false otherwise.
    bool less_than(const pos_t& a, const pos_t& b) const;

    /// Return True if position A is greater than position B in our sort, and false otherwise.
    bool greater_than(const Position& a, const Position& b) const;

  private:
    /// What's the maximum size of reads in serialized, uncompressed bytes to
    /// load into memory for a single temp file chunk, during the streaming
    /// sort?
    /// For reference, a whole-genome GAM file is about 500 GB of uncompressed data
    size_t max_buf_size = (512 * 1024 * 1024);
    /// What's the max fan-in when combining temp files, during the streaming sort?
    /// This will be computed based on the max file descriptor limit from the OS.
    size_t max_fan_in;
    
    using cursor_t = stream::ProtobufIterator<Alignment>;
    using emitter_t = stream::ProtobufEmitter<Alignment>;
    
    /// Open all the given input files, keeping the streams and cursors in the given lists.
    /// We use lists because none of these should be allowed to move after creation.
    void open_all(const vector<string>& filenames, list<ifstream>& streams, list<cursor_t>& cursors);
    
    /// Merge all the reads from the given list of cursors into the given emitter.
    /// The total expected number of reads can be passed for progress bar purposes.
    void streaming_merge(list<cursor_t>& cursors, emitter_t& emitter, size_t expected_reads = 0);
    
    /// Merge all the given temp input files into one or more temp output
    /// files, opening no more than max_fan_in input files at a time. The input
    /// files, which must be from temp_file::create(), will be deleted.
    ///
    /// If reads_per_file is specified, it will be used to show progress bars,
    /// and will be updated for newly-created files.
    vector<string> streaming_merge(const vector<string>& temp_names_in, unordered_map<string, size_t>* reads_per_file = nullptr);
};
}
#endif

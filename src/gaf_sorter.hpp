#ifndef VG_GAF_SORTER_HPP_INCLUDED
#define VG_GAF_SORTER_HPP_INCLUDED

/** \file 
 * Tools for sorting GAF records.
 *
 * NOTE: There is an equivalent standalone tool included with GBZ-base.
 *
 * TODO: Asynchronous I/O.
 * TODO: Option for automatic detection of merge width to guarantee <= 2 rounds.
 * TODO: Option for giving approximate batch size in bytes.
 */

#include <algorithm>
#include <cstdint>
#include <functional>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

#include <gbwt/utils.h>

namespace vg {

//------------------------------------------------------------------------------

/**
 * A record corresponding to a single line (alignment) in a GAF file.
 * The record contains an integer key and the original line.
 * Various types of keys can be derived from the value, but the line is not
 * parsed beyond that.
 */
struct GAFSorterRecord {
    /// Integer key.
    std::uint64_t key;

    /// GAF line.
    std::string value;

    /// Hasher used for random shuffling.
    static std::hash<std::string> hasher;

    /// Missing key. Records without a key are sorted to the end.
    constexpr static std::uint64_t MISSING_KEY = std::numeric_limits<std::uint64_t>::max();

    /// Types of keys that can be derived from the value.
    enum key_type {
        /// (minimum node id, maximum node id) in the path.
        key_node_interval,
        /// Hash of the value for random shuffling.
        key_hash,
    };

    /// Default constructor.
    GAFSorterRecord() : key(MISSING_KEY) {}

    /// Constructor that consumes the given value and sets the key.
    GAFSorterRecord(std::string&& value, key_type type) : key(MISSING_KEY), value(std::move(value)) {
        this->set_key(type);
    }

    /// Records are sorted by key in ascending order.
    bool operator<(const GAFSorterRecord& another) const {
        return (this->key < another.key);
    }

    /// Flips they key to reverse the order.
    /// Sorting is based on ascending order, while priority queues return the largest element first.
    void flip_key() {
        this->key = std::numeric_limits<std::uint64_t>::max() - this->key;
    }

    /// Sets a key of the given type, or MISSING_KEY if the key cannot be derived.
    void set_key(key_type type);

    /// Serializes the record to a stream. Returns true on success.
    bool serialize(std::ostream& out) const;

    /// Writes the underlying GAF line to a stream. Returns true on success.
    bool write_line(std::ostream& out) const;

    /// Deserializes the record from a stream. Returns true on success.
    bool deserialize(std::istream& in);

    /// Reads a GAF line from a stream and sets the key. Returns true on success.
    bool read_line(std::istream& in, key_type type);

    /// Returns a view of the given 0-based field, or an empty string if the field is missing.
    std::string_view get_field(size_t field) const;

    /// Calls the given function with a 0-based field index and the field value.
    /// Stops if the function returns false.
    void for_each_field(const std::function<bool(size_t, std::string_view)>& lambda) const;

    /// Returns the path as a GBWT path in forward orientation.
    /// If an ok flag is given, sets it to false and prints an error message
    /// if the path could not be parsed.
    gbwt::vector_type as_gbwt_path(bool* ok = nullptr) const;

    /// 0-based field number for the strand/orientation in a GAF line.
    constexpr static size_t STRAND_FIELD = 4;

    /// 0-based field number for the path in a GAF line.
    constexpr static size_t PATH_FIELD = 5;

    /// Number of mandatory fields in a GAF line.
    constexpr static size_t MANDATORY_FIELDS = 12;
};

//------------------------------------------------------------------------------

/**
 * A file of GAFSorterRecords or GAF lines.
 *
 * The records are sorted in increasing order by key.
 * The object is movable but not copyable.
 */
struct alignas(128) GAFSorterFile {
    /// File name.
    std::string name;

    /// Header lines that should be written to a raw GAF file.
    std::unique_ptr<std::vector<std::string>> header_lines;

    /// File name for the GBWT index, if any.
    std::string gbwt_file;

    /// Should the GBWT index be bidirectional?
    bool bidirectional_gbwt;

    /// Number of records.
    size_t records;

    /// Is this a temporary file created with temp_file::create()?
    bool temporary;

    /// Is this file compressed?
    bool compressed;

    /// Is this a raw GAF file?
    bool raw_gaf;

    /// Has the file been removed?
    bool removed;

    /// Success flag.
    bool ok;

    /// Default constructor that creates a compressed temporary file.
    GAFSorterFile();

    /// Constructor that creates a raw GAF file with the given name, header, and a GBWT index that may be bidirectional.
    /// Takes ownership of the header lines.
    GAFSorterFile(
        const std::string& name,
        std::unique_ptr<std::vector<std::string>> header_lines,
        const std::string& gbwt_file, bool bidirectional_gbwt
    );

    /// If the file is temporary, the destructor removes the file.
    ~GAFSorterFile();

    GAFSorterFile(const GAFSorterFile&) = delete;
    GAFSorterFile& operator=(const GAFSorterFile&) = delete;
    GAFSorterFile(GAFSorterFile&&) = default;
    GAFSorterFile& operator=(GAFSorterFile&&) = default;

    /// Returns an output stream to the file.
    /// The first return value is the actual stream.
    /// The second return value is a unique pointer which may contain a newly created stream.
    /// If this is a raw GAF file, header lines are written to the output.
    /// Sets the success flag.
    std::pair<std::ostream*, std::unique_ptr<std::ostream>> open_output();

    /// Writes the record to the file.
    /// Updates the number of records and sets the success flag.
    void write(const GAFSorterRecord& record, std::ostream& out) {
        this->ok &= (this->raw_gaf ? record.write_line(out) : record.serialize(out));
        this->records++;
    }

    /// Returns an input stream to the file.
    /// The first return value is the actual stream.
    /// The second return value is a unique pointer which may contain a newly created stream.
    /// Sets the success flag.
    std::pair<std::istream*, std::unique_ptr<std::istream>> open_input();

    /// Reads the next record from the file, assuming that this is not a raw GAF file.
    /// Sets the success flag.
    void read(GAFSorterRecord& record, std::istream& in) {
        this->ok &= (this->raw_gaf ? false : record.deserialize(in));
    }

    /// Returns true if the file is actually stdin/stdout.
    /// In that case, open_input() should not be called.
    bool is_std_in_out() const {
        return (this->name == "-");
    }

    /// Removes the file if it is temporary.
    void remove_temporary();
};

//------------------------------------------------------------------------------

/**
 * Parameters for the GAF sorter.
 */
struct GAFSorterParameters {
    /// Default for threads.
    constexpr static size_t THREADS = 1;

    /// Default for records_per_file.
    constexpr static size_t RECORDS_PER_FILE = 1000000;

    /// Default for files_per_merge.
    constexpr static size_t FILES_PER_MERGE = 32;

    /// Default for buffer_size.
    constexpr static size_t BUFFER_SIZE = 1000;

    /// Key type used for sorting.
    GAFSorterRecord::key_type key_type = GAFSorterRecord::key_node_interval;

    /// Number of parallel sort/merge jobs.
    size_t threads = THREADS;

    /// Number of records per file in the initial sort.
    size_t records_per_file = RECORDS_PER_FILE;

    /// Number of files to merge at once.
    size_t files_per_merge = FILES_PER_MERGE;

    /// Buffer size for reading and writing records.
    size_t buffer_size = BUFFER_SIZE;

    /// Use stable sorting.
    bool stable = false;

    /// GBWT output file, if any.
    std::string gbwt_file;

    /// Make the GBWT bidirectional.
    bool bidirectional_gbwt = false;

    /// Print progress information to stderr.
    bool progress = false;
};

/**
 * Sorts the given GAF file into the given output file.
 *
 * The initial round sorts the records into temporary files with params.records_per_file records each.
 * Each successive round merges the temporary files into larger files until there is only one file left.
 * Each merge job merges params.files_per_merge files.
 * Use "-" for reading stdin / writing to stdout.
 * Returns false and prints an error message on failure.
 *
 * If a GBWT file is specified in the parameters, a GBWT index of the paths is built for the paths and written to the file.
 * The GBWT index can be unidirectional or bidirectional, depending on the parameters.
 * Paths in the GBWT are in the same order as the records in the output file.
 * The GBWT index does not contain any metadata.
 * GBWT construction uses one additional thread.
 */
bool sort_gaf(const std::string& input_file, const std::string& output_file, const GAFSorterParameters& params);

/**
 * Sorts the given GAF lines into the given output file, with an option to use stable sorting.
 *
 * The lines are converted into GAFSorterRecord objects, with the given key type.
 * The original lines are consumed.
 * Sets the ok flag in the output and prints an error message on failure.
 *
 * If a GBWT file is specified in the output file, a GBWT index of the paths is built for the paths and written to the file.
 * The GBWT index can be unidirectional or bidirectional, depending on the output file.
 * Paths in the GBWT are in the same order as the records in the output file.
 * The GBWT index does not contain any metadata.
 * GBWT construction uses one additional thread.
 *
 * This function is intended to be used with std::thread.
 */
void sort_gaf_lines(std::unique_ptr<std::vector<std::string>> lines, GAFSorterRecord::key_type key_type, bool stable, GAFSorterFile& output);

/**
 * Merges the given files into a single output file.
 *
 * The records in each input file are assumed to be sorted with sort_gaf_lines().
 * Records are read and written in blocks of the given size.
 * If the input files are in the same order as the corresponding batches in the initial sort, this is a stable merge.
 * Consumes the inputs and removes the files if they are temporary.
 * Sets the ok flag in the output and prints an error message on failure.
 *
 * If a GBWT file is specified in the output file, a GBWT index of the paths is built for the paths and written to the file.
 * The GBWT index can be unidirectional or bidirectional, depending on the output file.
 * Paths in the GBWT are in the same order as the records in the output file.
 * The GBWT index does not contain any metadata.
 * GBWT construction uses one additional thread.
 *
 * This function is intended to be used with std::thread.
 */
void merge_gaf_records(std::unique_ptr<std::vector<GAFSorterFile>> inputs, GAFSorterFile& output, size_t buffer_size);

//------------------------------------------------------------------------------

} // namespace vg

#endif // VG_GAF_SORTER_HPP_INCLUDED

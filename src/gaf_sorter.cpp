#include "gaf_sorter.hpp"

#include <chrono>
#include <deque>
#include <fstream>
#include <queue>
#include <thread>
#include <charconv>

// Needed for the temporary file creation.
#include "utility.hpp"

// For reading and writing compressed temporary files.
#include "zstdutil.hpp"

// For reading compressed input.
#include <htslib/hfile.h>
#include <htslib/hts.h>

// For building a GBWT index of the paths.
#include <gbwt/dynamic_gbwt.h>

namespace vg {

//------------------------------------------------------------------------------

// Public class constants.

constexpr std::uint64_t GAFSorterRecord::MISSING_KEY;

constexpr size_t GAFSorterRecord::STRAND_FIELD;
constexpr size_t GAFSorterRecord::PATH_FIELD;
constexpr size_t GAFSorterRecord::MANDATORY_FIELDS;

constexpr size_t GAFSorterParameters::THREADS;
constexpr size_t GAFSorterParameters::RECORDS_PER_FILE;
constexpr size_t GAFSorterParameters::FILES_PER_MERGE;
constexpr size_t GAFSorterParameters::BUFFER_SIZE;

//------------------------------------------------------------------------------

void GAFSorterRecord::set_key(key_type type) {
    if (type == key_node_interval) {
        std::uint32_t min_id = std::numeric_limits<std::uint32_t>::max();
        std::uint32_t max_id = 0;
        std::string_view path = this->get_field(PATH_FIELD);
        size_t start = 1;
        while (start < path.size()) {
            std::uint32_t id = 0;
            auto result = std::from_chars(path.data() + start, path.data() + path.size(), id);
            if (result.ec != std::errc()) {
                this->key = MISSING_KEY;
                return;
            }
            min_id = std::min(min_id, id);
            max_id = std::max(max_id, id);
            start = (result.ptr - path.data()) + 1;
        }
        if (min_id == std::numeric_limits<std::uint32_t>::max()) {
            this->key = MISSING_KEY;
        } else {
            this->key = (static_cast<std::uint64_t>(min_id) << 32) | max_id;
        }
    } else if (type == key_hash) {
        this->key = hasher(this->value);
    } else {
        this->key = MISSING_KEY;
    }
}

bool GAFSorterRecord::serialize(std::ostream& out) const {
    bool success = true;
    out.write(reinterpret_cast<const char*>(&this->key), sizeof(this->key));
    success &= out.good();
    std::uint64_t length = this->value.size();
    out.write(reinterpret_cast<const char*>(&length), sizeof(length));
    success &= out.good();
    out.write(this->value.data(), length);
    success &= out.good();
    return success;
}

bool GAFSorterRecord::write_line(std::ostream& out) const {
    out << this->value << '\n';
    return out.good();
}

bool GAFSorterRecord::deserialize(std::istream& in) {
    bool success = true;
    in.read(reinterpret_cast<char*>(&this->key), sizeof(this->key));
    success &= in.good();
    std::uint64_t length = 0;
    in.read(reinterpret_cast<char*>(&length), sizeof(length));
    success &= in.good();
    this->value.resize(length);
    in.read(&this->value[0], length);
    success &= in.good();
    return success;
}

bool GAFSorterRecord::read_line(std::istream& in, key_type type) {
    std::getline(in, this->value);
    if (in.eof()) {
        return false;
    }
    this->set_key(type);
    return true;
}

std::string_view GAFSorterRecord::get_field(size_t field) const {
    std::string_view result;
    this->for_each_field([&](size_t i, std::string_view value) -> bool {
        if (i == field) {
            result = value;
            return false;
        }
        return true;
    });
    return result;
}

void GAFSorterRecord::for_each_field(const std::function<bool(size_t, std::string_view)>& lambda) const {
    size_t start = 0, end = 0;
    size_t i = 0;
    while (end != std::string::npos) {
        end = this->value.find('\t', start);
        if (!lambda(i, std::string_view(this->value).substr(start, end - start))) {
            break;
        }
        start = end + 1;
        i++;
    }
}

gbwt::vector_type GAFSorterRecord::as_gbwt_path(bool* ok) const {
    std::string_view strand, path;
    this->for_each_field([&](size_t i, std::string_view value) -> bool {
        if (i == STRAND_FIELD) {
            strand = value;
        } else if (i == PATH_FIELD) {
            path = value;
            return false; // Stop after the path.
        }
        return true; // Continue to the next field.
    });

    gbwt::vector_type result;
    if (path.size() == 1 && path[0] == '*') {
        // Unaligned read.
        return result;
    }

    size_t start = 0;
    while (start < path.size()) {
        bool is_reverse;
        if (path[start] == '<') {
            is_reverse = true;
        } else if (path[start] == '>') {
            is_reverse = false;
        } else {
            is_reverse = false;
            if (ok != nullptr) {
                *ok = false;
                std::cerr << "error: [gaf_sorter] invalid path: " << path << std::endl;
            }
            return result;
        }
        start++;
        gbwt::size_type node_id = 0;
        auto res = std::from_chars(path.data() + start, path.data() + path.size(), node_id);
        if (res.ec != std::errc() || node_id == 0) {
            if (ok != nullptr) {
                *ok = false;
                std::cerr << "error: [gaf_sorter] invalid path: " << path << std::endl;
            }
            return result;
        }
        result.push_back(gbwt::Node::encode(node_id, is_reverse));
        start = res.ptr - path.data();
    }

    // Now check the orientation.
    if (strand.size() == 1 && strand[0] == '-') {
        gbwt::reversePath(result);
    }

    return result;
}

//------------------------------------------------------------------------------

GAFSorterFile::GAFSorterFile() :
    name(temp_file::create("gaf-sorter")), header_lines(nullptr),
    gbwt_file(""), bidirectional_gbwt(false), records(0),
    temporary(true), compressed(true), raw_gaf(false), removed(false), ok(true) {
}

GAFSorterFile::GAFSorterFile(
    const std::string& name,
    std::unique_ptr<std::vector<std::string>> header_lines,
    const std::string& gbwt_file, bool bidirectional_gbwt
) :
    name(name), header_lines(std::move(header_lines)),
    gbwt_file(gbwt_file), bidirectional_gbwt(bidirectional_gbwt), records(0),
    temporary(false), compressed(false), raw_gaf(true), removed(false), ok(true) {
}

GAFSorterFile::~GAFSorterFile() {
    this->remove_temporary();
}

std::pair<std::ostream*, std::unique_ptr<std::ostream>> GAFSorterFile::open_output() {
    std::pair<std::ostream*, std::unique_ptr<std::ostream>> result;
    if (this->is_std_in_out()) {
        result.first = &std::cout;
    } else if (this->compressed) {
        result.second.reset(new zstd_ofstream(this->name));
        result.first = result.second.get();
    } else {
        result.second.reset(new std::ofstream(this->name, std::ios::binary));
        result.first = result.second.get();
    }
    if (!result.first->good()) {
        this->ok = false;
        std::cerr << "error: [gaf_sorter] could not open output file " << this->name << std::endl;
    }

    if (this->raw_gaf && this->header_lines != nullptr) {
        for (const std::string& line : *this->header_lines) {
            (*result.first) << line << '\n';
        }
        if (!result.first->good()) {
            this->ok = false;
            std::cerr << "error: [gaf_sorter] could not write header lines to output file " << this->name << std::endl;
        }
    }

    return result;
}

std::pair<std::istream*, std::unique_ptr<std::istream>> GAFSorterFile::open_input() {
    std::pair<std::istream*, std::unique_ptr<std::istream>> result;
    if (this->is_std_in_out()) {
        result.first = &std::cin;
    } else if (this->compressed) {
        result.second.reset(new zstd_ifstream(this->name));
        result.first = result.second.get();
    } else {
        result.second.reset(new std::ifstream(this->name, std::ios::binary));
        result.first = result.second.get();
    }
    if (!result.first->good()) {
        this->ok = false;
        std::cerr << "error: [gaf_sorter] could not open input file " << this->name << std::endl;
    }
    return result;
}

void GAFSorterFile::remove_temporary() {
    if (this->temporary && !this->removed) {
        temp_file::remove(this->name);
        this->removed = true;
        this->ok = false;
    }
}

//------------------------------------------------------------------------------

std::unique_ptr<gbwt::GBWTBuilder> create_gbwt_builder(const GAFSorterFile& output) {
    if (output.gbwt_file.empty()) {
        return nullptr;
    }
    gbwt::size_type node_width = sizeof(gbwt::vector_type::value_type) * 8;
    return std::make_unique<gbwt::GBWTBuilder>(node_width);
}

void finish_gbwt_construction(gbwt::GBWTBuilder* builder, GAFSorterFile& output) {
    if (builder == nullptr) {
        return;
    }
    builder->finish();
    try {
        sdsl::simple_sds::serialize_to(builder->index, output.gbwt_file);
    } catch (const std::runtime_error& e) {
        std::cerr << "error: [gaf_sorter] could not write GBWT index to " << output.gbwt_file << ": " << e.what() << std::endl;
        output.ok = false;
    }
}

//------------------------------------------------------------------------------

bool sort_gaf(const std::string& input_file, const std::string& output_file, const GAFSorterParameters& params) {
    // Timestamp for the start and the total number of records.
    auto start_time = std::chrono::high_resolution_clock::now();
    size_t total_records = 0;
    auto report_time = [&]() {
        if (params.progress) {
            auto end_time = std::chrono::high_resolution_clock::now();
            double seconds = std::chrono::duration<double>(end_time - start_time).count();
            std::cerr << "Sorted " << total_records << " records in " << seconds << " seconds" << std::endl;
        }
    };

    // Worker threads.
    size_t num_threads = std::max(params.threads, size_t(1));
    if (params.progress) {
        std::cerr << "Sorting GAF records with " << num_threads << " worker threads" << std::endl;
    }
    std::vector<std::thread> threads(num_threads);
    auto join_all = [&]() {
        for (std::thread& thread : threads) {
            if (thread.joinable()) {
                thread.join();
            }
        }
    };

    // Temporary output files.
    std::vector<std::unique_ptr<GAFSorterFile>> files;
    auto check_all_files = [&]() -> bool {
        bool all_ok = true;
        for (std::unique_ptr<GAFSorterFile>& file : files) {
            all_ok &= file->ok;
        }
        return all_ok;
    };

    // Open the input file and read header lines.
    std::unique_ptr<std::vector<std::string>> header_lines(new std::vector<std::string>());
    std::string peek;
    htsFile* input = hts_open(input_file.c_str(), "r");
    if (input == nullptr) {
        std::cerr << "error: [gaf_sorter] could not open input file " << input_file << std::endl;
        return false;
    }
    {
        kstring_t s_buffer = KS_INITIALIZE;
        while (hts_getline(input, '\n', &s_buffer) >= 0) {
            std::string line(ks_str(&s_buffer), ks_len(&s_buffer));
            if (!line.empty() && line[0] == '@') {
                header_lines->push_back(line);
            } else {
                peek = line;
                break;
            }
        }
        ks_free(&s_buffer);
    }

    // Initial sort. If a worker thread fails, we break on join.
    size_t batch = 0;
    size_t initial_batch_size = std::max(params.records_per_file, size_t(1));
    if (params.progress) {
        std::cerr << "Initial sort: " << initial_batch_size << " records per file" << std::endl;
    }
    while (true) {
        // Read the next batch.
        std::unique_ptr<std::vector<std::string>> lines(new std::vector<std::string>());
        lines->reserve(initial_batch_size);
        if (!peek.empty()) {
            lines->push_back(std::move(peek));
            peek.clear();
        }
        kstring_t s_buffer = KS_INITIALIZE;
        std::string line;
        while (lines->size() < initial_batch_size && hts_getline(input, '\n', &s_buffer) >= 0) {
            lines->push_back(std::string(ks_str(&s_buffer), ks_len(&s_buffer)));
        }
        total_records += lines->size();

        // Peek at the first line of the next batch to determine if there is only one batch.
        if (batch == 0) {
            if (hts_getline(input, '\n', &s_buffer) < 0) {
                if (params.progress) {
                    std::cerr << "Sorting directly to the final output" << std::endl;
                    if (!params.gbwt_file.empty()) {
                        std::cerr << "Building a GBWT index of the paths to " << params.gbwt_file << std::endl;
                    }
                }
                GAFSorterFile out(output_file, std::move(header_lines), params.gbwt_file, params.bidirectional_gbwt);
                sort_gaf_lines(std::move(lines), params.key_type, params.stable, out);
                ks_free(&s_buffer);
                hts_close(input);
                if (out.ok) {
                    report_time();
                    return true;
                } else {
                    return false;
                }
            }
            peek = std::string(ks_str(&s_buffer), ks_len(&s_buffer));
        }
        ks_free(&s_buffer);
        if (lines->empty()) {
            break;
        }

        // Sort the batch to a temporary file.
        std::unique_ptr<GAFSorterFile> out(new GAFSorterFile());
        size_t thread_id = batch % num_threads;
        if (threads[thread_id].joinable()) {
            threads[thread_id].join();
            if (!files[batch - num_threads]->ok) {
                break;
            }
        }
        threads[thread_id] = std::thread(sort_gaf_lines, std::move(lines), params.key_type, params.stable, std::ref(*out));
        files.push_back(std::move(out));
        batch++;
    }
    hts_close(input);
    join_all();
    if (!check_all_files()) {
        return false;
    }
    if (params.progress) {
        std::cerr << "Initial sort finished with " << total_records << " records in " << files.size() << " files" << std::endl;
    }

    // Intermediate merges. If a worker thread fails, we break on join.
    size_t files_per_merge = std::max(params.files_per_merge, size_t(2));
    size_t round = 0;
    while (files.size() > files_per_merge) {
        if (params.progress) {
            std::cerr << "Round " << round << ": " << files_per_merge << " files per batch" << std::endl;
        }
        std::vector<std::unique_ptr<GAFSorterFile>> next_files;
        batch = 0;
        for (size_t i = 0; i < files.size(); i += files_per_merge) {
            if (i + 1 == files.size()) {
                // If we have a single file left, just move it to the next round.
                next_files.push_back(std::move(files[i]));
                continue;
            }
            std::unique_ptr<std::vector<GAFSorterFile>> batch_files(new std::vector<GAFSorterFile>());
            for (size_t j = 0; j < files_per_merge && i + j < files.size(); j++) {
                batch_files->push_back(std::move(*(files[i + j])));
            }
            std::unique_ptr<GAFSorterFile> out(new GAFSorterFile());
            size_t thread_id = batch % num_threads;
            if (threads[thread_id].joinable()) {
                threads[thread_id].join();
                if (!next_files[batch - num_threads]->ok) {
                    break;
                }
            }
            threads[thread_id] = std::thread(merge_gaf_records, std::move(batch_files), std::ref(*out), params.buffer_size);
            next_files.push_back(std::move(out));
            batch++;
        }
        join_all();
        files = std::move(next_files);
        if (!check_all_files()) {
            return false;
        }
        if (params.progress) {
            std::cerr << "Round " << round << " finished with " << files.size() << " files" << std::endl;
        }
        round++;
    }

    // Final merge.
    {
        if (params.progress) {
            std::cerr << "Starting the final merge" << std::endl;
            if (!params.gbwt_file.empty()) {
                std::cerr << "Building a GBWT index of the paths to " << params.gbwt_file << std::endl;
            }
        }
        GAFSorterFile out(output_file, std::move(header_lines), params.gbwt_file, params.bidirectional_gbwt);
        std::unique_ptr<std::vector<GAFSorterFile>> inputs(new std::vector<GAFSorterFile>());
        for (std::unique_ptr<GAFSorterFile>& file : files) {
            inputs->push_back(std::move(*file));
        }
        merge_gaf_records(std::move(inputs), out, params.buffer_size);
        if (out.ok) {
            report_time();
            return true;
        } else {
            return false;
        }
    }
}

//------------------------------------------------------------------------------

void sort_gaf_lines(
    std::unique_ptr<std::vector<std::string>> lines,
    GAFSorterRecord::key_type key_type,
    bool stable,
    GAFSorterFile& output
) {
    if (lines == nullptr) {
        output.ok = false;
        std::cerr << "error: [gaf_sorter] sort_gaf_lines() called with null lines" << std::endl;
        return;
    }
    if (!output.ok || output.records > 0) {
        output.ok = false;
        std::cerr << "error: [gaf_sorter] sort_gaf_lines() called with an invalid output file" << std::endl;
        return;
    }

    // Convert the lines into GAFSorterRecord objects.
    std::vector<GAFSorterRecord> records;
    records.reserve(lines->size());
    for (std::string& line : *lines) {
        records.emplace_back(std::move(line), key_type);
    }
    lines.reset();

    // Sort the records.
    if (stable) {
        std::stable_sort(records.begin(), records.end());
    } else {
        std::sort(records.begin(), records.end());
    }

    // Create a GBWT index, if GBWT output is requested.
    std::unique_ptr<gbwt::GBWTBuilder> builder = create_gbwt_builder(output);

    // Write the sorted records to the output file.
    auto out = output.open_output();
    if (!output.ok) {
        // open_output() already prints an error message.
        return;
    }
    for (GAFSorterRecord& record : records) {
        output.write(record, *out.first);
        if (builder != nullptr) {
            gbwt::vector_type path = record.as_gbwt_path(&output.ok);
            if (!output.ok) {
                // We already printed an error message.
                return;
            }
            builder->insert(path, output.bidirectional_gbwt);
        }
    }
    out.first->flush(); // Just in case, if we are for example writing to std::cout.
    out.second.reset();

    // Finish the GBWT construction, if necessary.
    finish_gbwt_construction(builder.get(), output);
}

//------------------------------------------------------------------------------

void merge_gaf_records(std::unique_ptr<std::vector<GAFSorterFile>> inputs, GAFSorterFile& output, size_t buffer_size) {
    if (inputs == nullptr) {
        output.ok = false;
        std::cerr << "error: [gaf_sorter] merge_gaf_records() called with null inputs" << std::endl;
        return;
    }
    if (buffer_size == 0) {
        buffer_size = 1;
    }
    for (GAFSorterFile& input : *inputs) {
        if (!input.ok || input.raw_gaf) {
            output.ok = false;
            std::cerr << "error: [gaf_sorter] merge_gaf_records() called an invalid input file" << std::endl;
            return;
        }
    }
    if (!output.ok || output.records > 0) {
        output.ok = false;
        std::cerr << "error: [gaf_sorter] merge_gaf_records called() with an invalid output file" << std::endl;
        return;
    }

    // Open the input files.
    std::vector<std::pair<std::istream*, std::unique_ptr<std::istream>>> in; in.reserve(inputs->size());
    std::vector<size_t> remaining; remaining.reserve(inputs->size());
    for (GAFSorterFile& input : *inputs) {
        in.emplace_back(input.open_input());
        remaining.push_back(input.records);
        if (!input.ok) {
            // open_input() already prints an error message.
            output.ok = false;
            return;
        }
    }

    // Open the output file.
    auto out = output.open_output();
    if (!output.ok) {
        // open_output() already prints an error message.
        return;
    }

    // Input buffers.
    std::vector<std::deque<GAFSorterRecord>> records;
    records.resize(in.size());
    auto read_buffer = [&](size_t i) {
        size_t count = std::min(buffer_size, remaining[i]);
        if (count > 0) {
            records[i].clear();
            for (size_t j = 0; j < count; j++) {
                records[i].emplace_back();
                (*inputs)[i].read(records[i].back(), *(in[i].first));
                records[i].back().flip_key(); // Flip for the priority queue.
            }
            remaining[i] -= count;
            if (!(*inputs)[i].ok) {
                output.ok = false;
            }
        }
    };
    for (size_t i = 0; i < in.size(); i++) {
        read_buffer(i);
    }
    if (!output.ok) {
        std::cerr << "error: [gaf_sorter] merge_gaf_records() failed to read the initial buffers" << std::endl;
        return;
    }

    // Output buffers.
    std::vector<GAFSorterRecord> buffer;
    std::unique_ptr<gbwt::GBWTBuilder> builder = create_gbwt_builder(output);
    buffer.reserve(buffer_size);
    auto write_buffer = [&]() {
        for (GAFSorterRecord& record : buffer) {
            output.write(record, *out.first);
            if (builder != nullptr) {
                gbwt::vector_type path = record.as_gbwt_path(&output.ok);
                if (!output.ok) {
                    // We already printed an error message.
                    return;
                }
                builder->insert(path, output.bidirectional_gbwt);
            }
        }
        buffer.clear();
    };

    // Merge loop.
    std::priority_queue<std::pair<GAFSorterRecord, size_t>> queue;
    for (size_t i = 0; i < records.size(); i++) {
        if (!records[i].empty()) {
            // We already flipped the key for the max heap in `read_buffer()`.
            // We also need to flip the source index to get a stable order.
            size_t source_index = records.size() - 1 - i;
            queue.emplace(std::move(records[i].front()), source_index);
            records[i].pop_front();
        }
    }
    while (!queue.empty()) {
        GAFSorterRecord record = std::move(queue.top().first);
        record.flip_key(); // Restore the original key.
        size_t source = records.size() - 1 - queue.top().second; // Flip the source index back.
        queue.pop();
        buffer.push_back(std::move(record));
        if (buffer.size() >= buffer_size) {
            write_buffer();
            if (!output.ok) {
                std::cerr << "error: [gaf_sorter] merge_gaf_records() failed to write the output buffer" << std::endl;
                return;
            }
        }
        if (records[source].empty()) {
            read_buffer(source);
            if (!output.ok) {
                std::cerr << "error: [gaf_sorter] merge_gaf_records() failed to read from " << (*inputs)[source].name << std::endl;
                return;
            }
        }
        if (!records[source].empty()) {
            // Max heap; we use flipped keys and source indexes.
            size_t source_index = records.size() - 1 - source;
            queue.emplace(std::move(records[source].front()), source_index);
            records[source].pop_front();
        }
    }
    if (!buffer.empty()) {
        write_buffer();
        if (!output.ok) {
            std::cerr << "error: [gaf_sorter] merge_gaf_records() failed to write the output buffer" << std::endl;
            return;
        }
    }

    // Close the files.
    for (size_t i = 0; i < in.size(); i++) {
        in[i].first = nullptr;
        in[i].second.reset();
    }
    out.first->flush(); // Just in case, if we are for example writing to std::cout.
    out.second.reset();
    finish_gbwt_construction(builder.get(), output);
}

//------------------------------------------------------------------------------

} // namespace vg

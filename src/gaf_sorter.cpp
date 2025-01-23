#include "gaf_sorter.hpp"

#include <algorithm>
#include <charconv>
#include <chrono>
#include <deque>
#include <fstream>
#include <queue>
#include <thread>

// Needed for the temporary file creation.
#include "utility.hpp"

namespace vg {

//------------------------------------------------------------------------------

// Public class constants.

constexpr std::uint64_t GAFSorterRecord::MISSING_KEY;
constexpr std::string_view GAFSorterRecord::GBWT_OFFSET_TAG;

//------------------------------------------------------------------------------

void GAFSorterRecord::set_key(key_type type) {
    if (type == key_node_interval) {
        std::uint32_t min_id = std::numeric_limits<std::uint32_t>::max();
        std::uint32_t max_id = 0;
        std::string_view path = this->get_field(PATH_FIELD);
        std::string_view::size_type start = 1;
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
    } else if (type == key_gbwt_pos) {
        std::uint32_t node_id = std::numeric_limits<std::uint32_t>::max();
        std::uint32_t offset = std::numeric_limits<std::uint32_t>::max();
        this->for_each_field([&](size_t i, std::string_view value) -> bool {
            if (i == PATH_FIELD && value.size() > 1) {
                auto result = std::from_chars(value.data() + 1, value.data() + value.size(), node_id);
                if (result.ec != std::errc()) {
                    return false;
                }
            } else if (i >= MANDATORY_FIELDS) {
                constexpr size_t TAG_SIZE = GBWT_OFFSET_TAG.size();
                if (value.size() > TAG_SIZE && value.substr(0, TAG_SIZE) == GBWT_OFFSET_TAG) {
                    auto result = std::from_chars(value.data() + TAG_SIZE, value.data() + value.size(), offset);
                    return false;
                }
            }
            return true;
        });
        if (node_id == std::numeric_limits<std::uint32_t>::max() || offset == std::numeric_limits<std::uint32_t>::max()) {
            // We either did not find both fields or failed to parse them.
            this->key = MISSING_KEY;
        } else {
            this->key = (static_cast<std::uint64_t>(node_id) << 32) | offset;
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
    std::string_view::size_type start = 0;
    std::string_view::size_type end = 0;
    size_t i = 0;
    while (end != std::string_view::npos) {
        end = this->value.find('\t', start);
        if (!lambda(i, std::string_view(this->value).substr(start, end - start))) {
            break;
        }
        start = end + 1;
        i++;
    }
}

//------------------------------------------------------------------------------

GAFSorterFile::GAFSorterFile() :
    records(0),
    temporary(true), compressed(true), raw_gaf(false), removed(false), ok(true) {
    this->name = temp_file::create("gaf-sorter");
}

GAFSorterFile::GAFSorterFile(const std::string& name) :
    name(name), records(0),
    temporary(false), compressed(false), raw_gaf(true), removed(false), ok(true) {
}

GAFSorterFile::~GAFSorterFile() {
    this->remove_temporary();
}

std::pair<std::ostream*, std::unique_ptr<std::ostream>> GAFSorterFile::open_output() {
    std::pair<std::ostream*, std::unique_ptr<std::ostream>> result;
    if (this->is_stdout()) {
        result.first = &std::cout;
    } else {
        result.second.reset(new std::ofstream(this->name, std::ios::binary));
        result.first = result.second.get();
    }
    this->ok = result.first->good();
    return result;
}

std::unique_ptr<std::istream> GAFSorterFile::open_input() {
    std::unique_ptr<std::istream> result(new std::ifstream(this->name, std::ios::binary));
    this->ok = result->good();
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

void sort_gaf(std::istream& input, const std::string& output_file, const GAFSorterParameters& params) {
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

    size_t num_threads = std::max(params.threads, size_t(1));
    if (params.progress) {
        std::cerr << "Sorting GAF records with " << num_threads << " worker threads" << std::endl;
    }
    std::vector<std::thread> threads(num_threads);

    // Initial sort.
    size_t batch = 0;
    std::vector<std::unique_ptr<GAFSorterFile>> files;
    size_t initial_batch_size = std::max(params.records_per_file, size_t(1));
    if (params.progress) {
        std::cerr << "Initial sort: " << initial_batch_size << " records per file" << std::endl;
    }
    std::string peek;
    while (input) {
        // Read the next batch.
        std::unique_ptr<std::vector<std::string>> lines(new std::vector<std::string>());
        lines->reserve(initial_batch_size);
        if (!peek.empty()) {
            lines->push_back(std::move(peek));
            peek.clear();
        }
        std::string line;
        while (lines->size() < initial_batch_size && std::getline(input, line)) {
            lines->push_back(std::move(line));
        }
        total_records += lines->size();

        // Peek at the first line of the next batch to determine if there is only one batch.
        if (batch == 0) {
            std::getline(input, peek);
            if (!input) {
                if (params.progress) {
                    std::cerr << "Sorting directly to the final output" << std::endl;
                }
                std::unique_ptr<GAFSorterFile> out(new GAFSorterFile(output_file));
                sort_gaf_lines(std::move(lines), params.key_type, params.stable, std::ref(*out));
                report_time();
                return;
            }
        }
        if (lines->empty()) {
            break;
        }

        // Sort the batch to a temporary file.
        std::unique_ptr<GAFSorterFile> out(new GAFSorterFile());
        size_t thread_id = batch % num_threads;
        if (threads[thread_id].joinable()) {
            threads[thread_id].join();
        }
        threads[thread_id] = std::thread(sort_gaf_lines, std::move(lines), params.key_type, params.stable, std::ref(*out));
        files.push_back(std::move(out));
        batch++;
    }
    for (size_t i = 0; i < num_threads; i++) {
        if (threads[i].joinable()) {
            threads[i].join();
        }
    }
    if (params.progress) {
        std::cerr << "Initial sort finished with " << total_records << " records in " << files.size() << " files" << std::endl;
    }

    // Intermediate merges.
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
            }
            threads[thread_id] = std::thread(merge_gaf_records, std::move(batch_files), std::ref(*out), params.buffer_size);
            next_files.push_back(std::move(out));
            batch++;
        }
        for (size_t i = 0; i < num_threads; i++) {
            if (threads[i].joinable()) {
                threads[i].join();
            }
        }
        files = std::move(next_files);
        if (params.progress) {
            std::cerr << "Round " << round << " finished with " << files.size() << " files" << std::endl;
        }
        round++;
    }

    // Final merge.
    {
        if (params.progress) {
            std::cerr << "Starting the final merge" << std::endl;
        }
        GAFSorterFile out(output_file);
        std::unique_ptr<std::vector<GAFSorterFile>> inputs(new std::vector<GAFSorterFile>());
        for (std::unique_ptr<GAFSorterFile>& file : files) {
            inputs->push_back(std::move(*file));
        }
        merge_gaf_records(std::move(inputs), out, params.buffer_size);
        report_time();
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
        return;
    }
    if (!output.ok || output.records > 0) {
        output.ok = false;
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

    // Write the sorted records to the output file.
    auto out = output.open_output();
    if (!output.ok) {
        return;
    }
    for (GAFSorterRecord& record : records) {
        output.write(record, *out.first);
    }
    out.second.reset();
}

//------------------------------------------------------------------------------

void merge_gaf_records(std::unique_ptr<std::vector<GAFSorterFile>> inputs, GAFSorterFile& output, size_t buffer_size) {
    if (inputs == nullptr) {
        output.ok = false;
        return;
    }
    if (buffer_size == 0) {
        buffer_size = 1;
    }
    for (GAFSorterFile& input : *inputs) {
        if (!input.ok || input.raw_gaf) {
            output.ok = false;
            return;
        }
    }
    if (!output.ok || output.records > 0) {
        output.ok = false;
        return;
    }

    // Open the input files.
    std::vector<std::unique_ptr<std::istream>> in; in.reserve(inputs->size());
    std::vector<size_t> remaining; remaining.reserve(inputs->size());
    for (GAFSorterFile& input : *inputs) {
        in.emplace_back(input.open_input());
        remaining.push_back(input.records);
        if (!input.ok) {
            output.ok = false;
            return;
        }
    }

    // Open the output file.
    auto out = output.open_output();
    if (!output.ok) {
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
                (*inputs)[i].read(records[i].back(), *(in[i]));
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
        return;
    }

    // Output buffer.
    std::vector<GAFSorterRecord> buffer;
    buffer.reserve(buffer_size);
    auto write_buffer = [&]() {
        for (GAFSorterRecord& record : buffer) {
            output.write(record, *out.first);
        }
        buffer.clear();
    };

    // Merge loop.
    std::priority_queue<std::pair<GAFSorterRecord, size_t>> queue;
    for (size_t i = 0; i < records.size(); i++) {
        if (!records[i].empty()) {
            queue.emplace(records[i].front(), i);
            records[i].pop_front();
        }
    }
    while (!queue.empty()) {
        GAFSorterRecord record = std::move(queue.top().first);
        record.flip_key(); // Restore the original key.
        size_t source = queue.top().second;
        queue.pop();
        buffer.push_back(std::move(record));
        if (buffer.size() >= buffer_size) {
            write_buffer();
            if (!output.ok) {
                return;
            }
        }
        if (records[source].empty()) {
            read_buffer(source);
            if (!output.ok) {
                return;
            }
        }
        if (!records[source].empty()) {
            queue.emplace(records[source].front(), source);
            records[source].pop_front();
        }
    }
    if (!buffer.empty()) {
        write_buffer();
    }

    // Close the files.
    for (size_t i = 0; i < in.size(); i++) {
        in[i].reset();
    }
    out.second.reset();
}

//------------------------------------------------------------------------------

} // namespace vg

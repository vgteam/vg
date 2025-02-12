/** \file
 *
 * Unit tests for gaf_sorter.cpp, which provides tools for sorting GAF records.
 */

#include "../gaf_sorter.hpp"
#include "../utility.hpp"

#include <fstream>
#include <random>

#include "catch.hpp"

namespace vg {

namespace unittest {

//------------------------------------------------------------------------------

namespace {

struct GAFInfo {
    std::string line;
    size_t id;
    std::uint32_t min_node, max_node;
    std::uint32_t first_node, gbwt_offset;
    bool first_is_reverse;

    GAFInfo(size_t id, size_t nodes, bool with_gbwt_offset) :
        id(id),
        min_node(std::numeric_limits<std::uint32_t>::max()), max_node(0),
        first_node(std::numeric_limits<std::uint32_t>::max()), gbwt_offset(std::numeric_limits<std::uint32_t>::max()),
        first_is_reverse(false) {

        // Name, query length, query start, query end, strand;
        std::string nd = std::to_string(nodes);
        this->line = "read" + std::to_string(id) + "\t" + nd + "\t0\t" + nd + "\t+\t";
 
        // Path.
        std::mt19937 rng(id);
        for (size_t i = 0; i < nodes; i++) {
            std::uint32_t node = rng() % 1000 + 1;
            bool reverse = rng() % 2;
            this->min_node = std::min(this->min_node, node);
            this->max_node = std::max(this->max_node, node);
            if (i == 0) {
                this->first_node = node;
                if (with_gbwt_offset) {
                    this->gbwt_offset = rng() % 1000;
                }
                this->first_is_reverse = reverse;
            }
            this->line += (reverse ? "<" : ">") + std::to_string(node);
        }

        // Path length, path start, path end, matches, alignment length, mapping quality.
        this->line += "\t" + nd + "\t0\t" + nd + "\t" + nd + "\t" + nd + "\t60";

        // Some arbitrary tags.
        this->line += "\tab:Z:cd\tef:i:42";

        // GBWT offset.
        if (with_gbwt_offset) {
            this->line.push_back('\t');
            this->line += GAFSorterRecord::GBWT_OFFSET_TAG;
            this->line += std::to_string(this->gbwt_offset);
        }

        // More tags.
        this->line += "\tgh:Z:ij\tkl:i:42";
    }

    std::uint64_t key(GAFSorterRecord::key_type type) const {
        if (type == GAFSorterRecord::key_node_interval) {
            return (static_cast<std::uint64_t>(this->min_node) << 32) | this->max_node;
        } else if (type == GAFSorterRecord::key_gbwt_pos) {
            if (this->gbwt_offset == std::numeric_limits<std::uint32_t>::max()) {
                return GAFSorterRecord::MISSING_KEY;
            } else {
                return (static_cast<std::uint64_t>(this->first_node) << 33) |
                    (static_cast<std::uint64_t>(this->first_is_reverse) << 32) | this->gbwt_offset;
            }
        } else if (type == GAFSorterRecord::key_hash) {
            return GAFSorterRecord::hasher(this->line);        
        } else {
            return GAFSorterRecord::MISSING_KEY;
        }
    }

    // Returns a copy of the line that can be consumed.
    std::string value() const {
        return this->line;
    }

    // Returns the id encoded in read name.
    static std::uint32_t decode_id(const std::string& line) {
        size_t pos = line.find('\t');
        return std::stoul(line.substr(4, pos - 4));
    }
};

std::unique_ptr<std::vector<std::string>> generate_gaf(size_t count, size_t path_length, double unaligned_probability) {
    std::unique_ptr<std::vector<std::string>> result(new std::vector<std::string>());
    result->reserve(count);
    std::mt19937 rng(count ^ path_length);
    for (size_t id = 0; id < count; id++) {
        double p = static_cast<double>(rng()) / rng.max();
        if (p < unaligned_probability) {
            GAFInfo info(id, 0, false);
            result->push_back(info.value());
        } else {
            GAFInfo info(id, path_length, true);
            result->push_back(info.value());
        }
    }
    return result;
}

std::vector<GAFSorterRecord> generate_records(size_t count, size_t path_length, double unaligned_probability) {
    auto lines = generate_gaf(count, path_length, unaligned_probability);
    std::vector<GAFSorterRecord> result;
    for (std::string line : *lines) {
        result.emplace_back(std::move(line), GAFSorterRecord::key_node_interval);
    }
    return result;
}

GAFSorterFile generate_sorted(size_t count, size_t path_length, double unaligned_probability, bool stable, const std::string* filename = nullptr) {
    auto lines = generate_gaf(count, path_length, unaligned_probability);
    GAFSorterFile output = (filename == nullptr ? GAFSorterFile() : GAFSorterFile(*filename));
    sort_gaf_lines(std::move(lines), GAFSorterRecord::key_node_interval, stable, output);
    return output;
}

void check_sorted(GAFSorterFile& file, bool raw_gaf, size_t lines, GAFSorterRecord::key_type key_type, bool stable) {
    REQUIRE(file.ok);
    REQUIRE(file.records == lines);

    std::pair<std::istream*, std::unique_ptr<std::istream>> in;
    if (raw_gaf) {
        in.second = std::unique_ptr<std::istream>(new std::ifstream(file.name, std::ios::binary));
        in.first = in.second.get();
    } else {
        in = file.open_input();
    }

    size_t line_num = 0;
    GAFSorterRecord previous;
    while (line_num < lines) {
        GAFSorterRecord record;
        if (raw_gaf) {
            REQUIRE(record.read_line(*in.first, key_type));
        } else {
            REQUIRE(record.deserialize(*in.first));
        }
        REQUIRE((line_num == 0 || previous.key <= record.key));
        if (stable && line_num > 0 && previous.key == record.key) {
            std::uint32_t prev_id = GAFInfo::decode_id(previous.value);
            std::uint32_t curr_id = GAFInfo::decode_id(record.value);
            REQUIRE(prev_id < curr_id);
        }
        previous = record;
        line_num++;
    }

    // There should not be any additional data.
    char c;
    REQUIRE(!in.first->get(c));
}

void merge_and_check(std::unique_ptr<std::vector<GAFSorterFile>> inputs, size_t buffer_size, size_t expected_records, GAFSorterRecord::key_type key_type) {
    std::string filename = temp_file::create("gaf-sorter");
    GAFSorterFile output(filename);
    merge_gaf_records(std::move(inputs), output, buffer_size);
    check_sorted(output, true, expected_records, key_type, false);
    temp_file::remove(filename);
}

void integrated_test(size_t count, size_t path_length, double unaligned_probability, const GAFSorterParameters& params) {
    // Generate the input.
    std::string input_file = temp_file::create("gaf-sorter");
    std::ofstream out(input_file, std::ios::binary);
    auto lines = generate_gaf(count, path_length, unaligned_probability);
    for (const std::string& line : *lines) {
        out << line << '\n';
    }
    out.close();

    // Sort the input.
    std::string output_file = temp_file::create("gaf-sorter");
    REQUIRE(sort_gaf(input_file, output_file, params));
    temp_file::remove(input_file);

    // Check the output.
    GAFSorterFile output(output_file);
    output.records = count; // This is a new file object, so we need to set the record count.
    check_sorted(output, true, count, params.key_type, false);
    temp_file::remove(output_file);
}

} // anonymous namespace

//------------------------------------------------------------------------------

TEST_CASE("Records and keys", "[gaf_sorter]") {
    SECTION("node interval keys") {
        for (size_t id = 0; id < 10; id++) {
            GAFInfo info(id, 10, false);
            GAFSorterRecord record(info.value(), GAFSorterRecord::key_node_interval);
            REQUIRE(record.key == info.key(GAFSorterRecord::key_node_interval));
        }
    }

    SECTION("node interval key with an empty path") {
        for (size_t id = 0; id < 10; id++) {
            GAFInfo info(id, 0, false);
            GAFSorterRecord record(info.value(), GAFSorterRecord::key_node_interval);
            REQUIRE(record.key == GAFSorterRecord::MISSING_KEY);
        }
    }

    SECTION("GBWT position keys") {
        for (size_t id = 0; id < 10; id++) {
            GAFInfo info(id, 10, true);
            GAFSorterRecord record(info.value(), GAFSorterRecord::key_gbwt_pos);
            REQUIRE(record.key == info.key(GAFSorterRecord::key_gbwt_pos));
        }
    }

    SECTION("GBWT position key with an empty path") {
        for (size_t id = 0; id < 10; id++) {
            GAFInfo info(id, 0, true);
            GAFSorterRecord record(info.value(), GAFSorterRecord::key_gbwt_pos);
            REQUIRE(record.key == GAFSorterRecord::MISSING_KEY);
        }
    }

    SECTION("GBWT position key without offset") {
        for (size_t id = 0; id < 10; id++) {
            GAFInfo info(id, 10, false);
            GAFSorterRecord record(info.value(), GAFSorterRecord::key_gbwt_pos);
            REQUIRE(record.key == GAFSorterRecord::MISSING_KEY);
        }
    }

    SECTION("hash keys") {
        for (size_t id = 0; id < 10; id++) {
            GAFInfo info(id, 10, false);
            GAFSorterRecord record(info.value(), GAFSorterRecord::key_hash);
            REQUIRE(record.key == info.key(GAFSorterRecord::key_hash));
        }
    }
}

TEST_CASE("Record serialization", "[gaf_sorter]") {
    SECTION("records") {
        auto records = generate_records(100, 10, 0.05);
        std::string filename = temp_file::create("gaf-sorter");
        std::ofstream out(filename, std::ios::binary);
        for (const GAFSorterRecord& record : records) {
            record.serialize(out);
        }
        out.close();

        std::ifstream in(filename, std::ios::binary);
        for (size_t i = 0; i < records.size(); i++) {
            GAFSorterRecord record;
            REQUIRE(record.deserialize(in));
            REQUIRE(record.key == records[i].key);
            REQUIRE(record.value == records[i].value);
        }
        char c;
        REQUIRE(!in.get(c));
    }

    SECTION("lines") {
        auto records = generate_records(120, 10, 0.05);
        std::string filename = temp_file::create("gaf-sorter");
        std::ofstream out(filename, std::ios::binary);
        for (const GAFSorterRecord& record : records) {
            record.write_line(out);
        }
        out.close();

        std::ifstream in(filename, std::ios::binary);
        for (size_t i = 0; i < records.size(); i++) {
            GAFSorterRecord record;
            REQUIRE(record.read_line(in, GAFSorterRecord::key_node_interval));
            REQUIRE(record.key == records[i].key);
            REQUIRE(record.value == records[i].value);
        }
        char c;
        REQUIRE(!in.get(c));
    }
}

//------------------------------------------------------------------------------

TEST_CASE("Sorting GAF records", "[gaf_sorter]") {
    SECTION("raw GAF output") {
        size_t n = 1000;
        std::string filename = temp_file::create("gaf-sorter");
        GAFSorterFile output = generate_sorted(n, 10, 0.05, false, &filename);
        check_sorted(output, true, n, GAFSorterRecord::key_node_interval, false);
        temp_file::remove(filename);
    }

    SECTION("empty GAF output") {
        size_t n = 0;
        std::string filename = temp_file::create("gaf-sorter");
        GAFSorterFile output = generate_sorted(n, 10, 0.05, false, &filename);
        check_sorted(output, true, n, GAFSorterRecord::key_node_interval, false);
        temp_file::remove(filename);
    }

    SECTION("record output") {
        size_t n = 1234;
        GAFSorterFile output = generate_sorted(n, 10, 0.05, false);
        check_sorted(output, false, n, GAFSorterRecord::key_node_interval, false);
    }

    SECTION("empty record output") {
        size_t n = 0;
        GAFSorterFile output = generate_sorted(n, 10, 0.05, false);
        check_sorted(output, false, n, GAFSorterRecord::key_node_interval, false);
    }

    SECTION("stable sorting") {
        size_t n = 1000;
        GAFSorterFile output = generate_sorted(n, 10, 0.05, true);
        check_sorted(output, false, n, GAFSorterRecord::key_node_interval, true);
    }
}

//------------------------------------------------------------------------------

TEST_CASE("Merging sorted files", "[gaf_sorter]") {
    SECTION("three files") {
        size_t n = 1000, expected_records = 0;
        std::unique_ptr<std::vector<GAFSorterFile>> inputs(new std::vector<GAFSorterFile>());
        for (size_t i = 0; i < 3; i++) {
            inputs->push_back(generate_sorted(n + i, 10, 0.05, false));
            expected_records += inputs->back().records;
        }
        merge_and_check(std::move(inputs), 100, expected_records, GAFSorterRecord::key_node_interval);
    }

    SECTION("one file is empty") {
        size_t n = 1000, expected_records = 0;
        std::unique_ptr<std::vector<GAFSorterFile>> inputs(new std::vector<GAFSorterFile>());
        for (size_t i = 0; i < 3; i++) {
            size_t count = (i == 1 ? 0 : n + i);
            inputs->push_back(generate_sorted(count, 10, 0.05, false));
            expected_records += inputs->back().records;
        }
        merge_and_check(std::move(inputs), 100, expected_records, GAFSorterRecord::key_node_interval);
    }

    SECTION("all files are empty") {
        size_t expected_records = 0;
        std::unique_ptr<std::vector<GAFSorterFile>> inputs(new std::vector<GAFSorterFile>());
        for (size_t i = 0; i < 3; i++) {
            inputs->push_back(generate_sorted(0, 10, 0.05, false));
            expected_records += inputs->back().records;
        }
        merge_and_check(std::move(inputs), 100, expected_records, GAFSorterRecord::key_node_interval);
    }

    SECTION("no input files") {
        size_t expected_records = 0;
        std::unique_ptr<std::vector<GAFSorterFile>> inputs(new std::vector<GAFSorterFile>());
        merge_and_check(std::move(inputs), 100, expected_records, GAFSorterRecord::key_node_interval);
    }
}

//------------------------------------------------------------------------------

TEST_CASE("GAF sorting", "[gaf_sorter]") {
    SECTION("one full batch") {
        size_t n = 1000;
        GAFSorterParameters params;
        params.records_per_file = 1000;
        integrated_test(n, 10, 0.05, params);
    }

    SECTION("one partial batch") {
        size_t n = 500;
        GAFSorterParameters params;
        params.records_per_file = 1000;
        integrated_test(n, 10, 0.05, params);
    }

    SECTION("one merge") {
        size_t n = 2000;
        GAFSorterParameters params;
        params.records_per_file = 1000;
        params.files_per_merge = 2;
        integrated_test(n, 10, 0.05, params);
    }

    SECTION("one merge + one batch") {
        size_t n = 3000;
        GAFSorterParameters params;
        params.records_per_file = 1000;
        params.files_per_merge = 2;
        integrated_test(n, 10, 0.05, params);
    }

    SECTION("multiple levels of merges") {
        size_t n = 10000;
        GAFSorterParameters params;
        params.records_per_file = 1000;
        params.files_per_merge = 2;
        integrated_test(n, 10, 0.05, params);
    }

    SECTION("multithreaded") {
        size_t n = 10000;
        GAFSorterParameters params;
        params.records_per_file = 1000;
        params.files_per_merge = 2;
        params.threads = 2;
        integrated_test(n, 10, 0.05, params);
    }

    SECTION("empty input") {
        size_t n = 0;
        GAFSorterParameters params;
        integrated_test(n, 10, 0.05, params);
    }
}

//------------------------------------------------------------------------------

} // namespace unittest

} // namespace vg

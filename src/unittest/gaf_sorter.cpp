/** \file
 *
 * Unit tests for gaf_sorter.cpp, which provides tools for sorting GAF records.
 */

#include "../gaf_sorter.hpp"
#include "../utility.hpp"

#include <fstream>
#include <random>

#include <gbwt/gbwt.h>

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

GAFSorterFile generate_sorted_temp(size_t count, size_t path_length, double unaligned_probability, bool stable) {
    auto lines = generate_gaf(count, path_length, unaligned_probability);
    GAFSorterFile output;
    sort_gaf_lines(std::move(lines), GAFSorterRecord::key_node_interval, stable, output);
    return output;
}

GAFSorterFile generate_sorted_gaf(size_t count, size_t path_length, double unaligned_probability, bool stable, const std::string& filename, const std::string& gbwt_file = "") {
    auto lines = generate_gaf(count, path_length, unaligned_probability);
    GAFSorterFile output(filename, gbwt_file);
    sort_gaf_lines(std::move(lines), GAFSorterRecord::key_node_interval, stable, output);
    return output;
}

gbwt::vector_type parse_path(const GAFSorterRecord& record) {
    gbwt::vector_type path;
    str_view field = record.get_field(GAFSorterRecord::PATH_FIELD);
    size_t start = 0;
    while (start < field.size) {
        bool is_reverse = (field[start] == '<');
        size_t end = start + 1;
        while (end < field.size && field[end] != '<' && field[end] != '>') {
            end++;
        }
        gbwt::size_type node_id = std::stoul(std::string(field.data + start + 1, field.data + end));
        path.push_back(gbwt::Node::encode(node_id, is_reverse));
        start = end;
    }
    return path;
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

    gbwt::GBWT gbwt_index;
    if (!file.gbwt_file.empty()) {
        sdsl::simple_sds::load_from(gbwt_index, file.gbwt_file);
        REQUIRE(gbwt_index.bidirectional());
        REQUIRE(gbwt_index.sequences() == 2 * lines);
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
        if (!file.gbwt_file.empty()) {
            gbwt::vector_type extracted = gbwt_index.extract(gbwt::Path::encode(line_num, false));
            gbwt::vector_type truth = parse_path(record);
            REQUIRE(extracted.size() == truth.size());
            for (size_t i = 0; i < extracted.size(); i++) {
                REQUIRE(extracted[i] == truth[i]);
            }
        }
        previous = record;
        line_num++;
    }

    // There should not be any additional data.
    char c;
    REQUIRE(!in.first->get(c));
}

void merge_and_check(std::unique_ptr<std::vector<GAFSorterFile>> inputs, size_t buffer_size, size_t expected_records, GAFSorterRecord::key_type key_type, bool with_gbwt) {
    std::string filename = temp_file::create("gaf-sorter");
    std::string gbwt_file = (with_gbwt ? temp_file::create("gaf-sorter") : "");
    GAFSorterFile output(filename, gbwt_file);
    merge_gaf_records(std::move(inputs), output, buffer_size);
    check_sorted(output, true, expected_records, key_type, false);
    temp_file::remove(filename);
    if (with_gbwt) {
        temp_file::remove(gbwt_file);
    }
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
    GAFSorterFile output(output_file, params.gbwt_file);
    output.records = count; // This is a new file object, so we need to set the record count.
    check_sorted(output, true, count, params.key_type, false);

    // Remove the output file(s).
    temp_file::remove(output_file);
    if (!params.gbwt_file.empty()) {
        temp_file::remove(params.gbwt_file);
    }
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
        GAFSorterFile output = generate_sorted_gaf(n, 10, 0.05, false, filename);
        check_sorted(output, true, n, GAFSorterRecord::key_node_interval, false);
        temp_file::remove(filename);
    }

    SECTION("raw GAF output with GBWT") {
        size_t n = 1000;
        std::string filename = temp_file::create("gaf-sorter");
        std::string gbwt_file = temp_file::create("gaf-sorter");
        GAFSorterFile output = generate_sorted_gaf(n, 10, 0.05, false, filename, gbwt_file);
        check_sorted(output, true, n, GAFSorterRecord::key_node_interval, false);
        temp_file::remove(filename);
        temp_file::remove(gbwt_file);
    }

    SECTION("empty GAF output") {
        size_t n = 0;
        std::string filename = temp_file::create("gaf-sorter");
        GAFSorterFile output = generate_sorted_gaf(n, 10, 0.05, false, filename);
        check_sorted(output, true, n, GAFSorterRecord::key_node_interval, false);
        temp_file::remove(filename);
    }

    SECTION("empty GAF output with GBWT") {
        size_t n = 0;
        std::string filename = temp_file::create("gaf-sorter");
        std::string gbwt_file = temp_file::create("gaf-sorter");
        GAFSorterFile output = generate_sorted_gaf(n, 10, 0.05, false, filename, gbwt_file);
        check_sorted(output, true, n, GAFSorterRecord::key_node_interval, false);
        temp_file::remove(filename);
        temp_file::remove(gbwt_file);
    }

    SECTION("record output") {
        size_t n = 1234;
        GAFSorterFile output = generate_sorted_temp(n, 10, 0.05, false);
        check_sorted(output, false, n, GAFSorterRecord::key_node_interval, false);
    }

    SECTION("empty record output") {
        size_t n = 0;
        GAFSorterFile output = generate_sorted_temp(n, 10, 0.05, false);
        check_sorted(output, false, n, GAFSorterRecord::key_node_interval, false);
    }

    SECTION("stable sorting") {
        size_t n = 1000;
        GAFSorterFile output = generate_sorted_temp(n, 10, 0.05, true);
        check_sorted(output, false, n, GAFSorterRecord::key_node_interval, true);
    }
}

//------------------------------------------------------------------------------

TEST_CASE("Merging sorted files", "[gaf_sorter]") {
    SECTION("three files") {
        size_t n = 1000, expected_records = 0;
        std::unique_ptr<std::vector<GAFSorterFile>> inputs(new std::vector<GAFSorterFile>());
        for (size_t i = 0; i < 3; i++) {
            inputs->push_back(generate_sorted_temp(n + i, 10, 0.05, false));
            expected_records += inputs->back().records;
        }
        merge_and_check(std::move(inputs), 100, expected_records, GAFSorterRecord::key_node_interval, false);
    }

    SECTION("three files with GBWT") {
        size_t n = 1000, expected_records = 0;
        std::unique_ptr<std::vector<GAFSorterFile>> inputs(new std::vector<GAFSorterFile>());
        for (size_t i = 0; i < 3; i++) {
            inputs->push_back(generate_sorted_temp(n + i, 10, 0.05, false));
            expected_records += inputs->back().records;
        }
        merge_and_check(std::move(inputs), 100, expected_records, GAFSorterRecord::key_node_interval, true);
    }

    SECTION("one file is empty") {
        size_t n = 1000, expected_records = 0;
        std::unique_ptr<std::vector<GAFSorterFile>> inputs(new std::vector<GAFSorterFile>());
        for (size_t i = 0; i < 3; i++) {
            size_t count = (i == 1 ? 0 : n + i);
            inputs->push_back(generate_sorted_temp(count, 10, 0.05, false));
            expected_records += inputs->back().records;
        }
        merge_and_check(std::move(inputs), 100, expected_records, GAFSorterRecord::key_node_interval, false);
    }

    SECTION("one file is empty with GBWT") {
        size_t n = 1000, expected_records = 0;
        std::unique_ptr<std::vector<GAFSorterFile>> inputs(new std::vector<GAFSorterFile>());
        for (size_t i = 0; i < 3; i++) {
            size_t count = (i == 1 ? 0 : n + i);
            inputs->push_back(generate_sorted_temp(count, 10, 0.05, false));
            expected_records += inputs->back().records;
        }
        merge_and_check(std::move(inputs), 100, expected_records, GAFSorterRecord::key_node_interval, true);
    }

    SECTION("all files are empty") {
        size_t expected_records = 0;
        std::unique_ptr<std::vector<GAFSorterFile>> inputs(new std::vector<GAFSorterFile>());
        for (size_t i = 0; i < 3; i++) {
            inputs->push_back(generate_sorted_temp(0, 10, 0.05, false));
            expected_records += inputs->back().records;
        }
        merge_and_check(std::move(inputs), 100, expected_records, GAFSorterRecord::key_node_interval, false);
    }

    SECTION("all files are empty with GBWT") {
        size_t expected_records = 0;
        std::unique_ptr<std::vector<GAFSorterFile>> inputs(new std::vector<GAFSorterFile>());
        for (size_t i = 0; i < 3; i++) {
            inputs->push_back(generate_sorted_temp(0, 10, 0.05, false));
            expected_records += inputs->back().records;
        }
        merge_and_check(std::move(inputs), 100, expected_records, GAFSorterRecord::key_node_interval, true);
    }

    SECTION("no input files") {
        size_t expected_records = 0;
        std::unique_ptr<std::vector<GAFSorterFile>> inputs(new std::vector<GAFSorterFile>());
        merge_and_check(std::move(inputs), 100, expected_records, GAFSorterRecord::key_node_interval, false);
    }

    SECTION("no input files with GBWT") {
        size_t expected_records = 0;
        std::unique_ptr<std::vector<GAFSorterFile>> inputs(new std::vector<GAFSorterFile>());
        merge_and_check(std::move(inputs), 100, expected_records, GAFSorterRecord::key_node_interval, true);
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

    SECTION("one full batch with GBWT") {
        size_t n = 1000;
        GAFSorterParameters params;
        params.records_per_file = 1000;
        params.gbwt_file = temp_file::create("gaf-sorter");
        integrated_test(n, 10, 0.05, params);
        temp_file::remove(params.gbwt_file);
    }

    SECTION("one partial batch") {
        size_t n = 500;
        GAFSorterParameters params;
        params.records_per_file = 1000;
        integrated_test(n, 10, 0.05, params);
    }

    SECTION("one partial batch with GBWT") {
        size_t n = 500;
        GAFSorterParameters params;
        params.records_per_file = 1000;
        params.gbwt_file = temp_file::create("gaf-sorter");
        integrated_test(n, 10, 0.05, params);
        temp_file::remove(params.gbwt_file);
    }

    SECTION("one merge") {
        size_t n = 2000;
        GAFSorterParameters params;
        params.records_per_file = 1000;
        params.files_per_merge = 2;
        integrated_test(n, 10, 0.05, params);
    }

    SECTION("one merge with GBWT") {
        size_t n = 2000;
        GAFSorterParameters params;
        params.records_per_file = 1000;
        params.files_per_merge = 2;
        params.gbwt_file = temp_file::create("gaf-sorter");
        integrated_test(n, 10, 0.05, params);
        temp_file::remove(params.gbwt_file);
    }

    SECTION("one merge + one batch") {
        size_t n = 3000;
        GAFSorterParameters params;
        params.records_per_file = 1000;
        params.files_per_merge = 2;
        integrated_test(n, 10, 0.05, params);
    }

    SECTION("one merge + one batch with GBWT") {
        size_t n = 3000;
        GAFSorterParameters params;
        params.records_per_file = 1000;
        params.files_per_merge = 2;
        params.gbwt_file = temp_file::create("gaf-sorter");
        integrated_test(n, 10, 0.05, params);
        temp_file::remove(params.gbwt_file);
    }

    SECTION("multiple levels of merges") {
        size_t n = 10000;
        GAFSorterParameters params;
        params.records_per_file = 1000;
        params.files_per_merge = 2;
        integrated_test(n, 10, 0.05, params);
    }

    SECTION("multiple levels of merges with GBWT") {
        size_t n = 10000;
        GAFSorterParameters params;
        params.records_per_file = 1000;
        params.files_per_merge = 2;
        params.gbwt_file = temp_file::create("gaf-sorter");
        integrated_test(n, 10, 0.05, params);
        temp_file::remove(params.gbwt_file);
    }

    SECTION("multithreaded") {
        size_t n = 10000;
        GAFSorterParameters params;
        params.records_per_file = 1000;
        params.files_per_merge = 2;
        params.threads = 2;
        integrated_test(n, 10, 0.05, params);
    }

    SECTION("multithreaded with GBWT") {
        size_t n = 10000;
        GAFSorterParameters params;
        params.records_per_file = 1000;
        params.files_per_merge = 2;
        params.threads = 2;
        params.gbwt_file = temp_file::create("gaf-sorter");
        integrated_test(n, 10, 0.05, params);
        temp_file::remove(params.gbwt_file);
    }

    SECTION("empty input") {
        size_t n = 0;
        GAFSorterParameters params;
        integrated_test(n, 10, 0.05, params);
    }

    SECTION("empty input with GBWT") {
        size_t n = 0;
        GAFSorterParameters params;
        params.gbwt_file = temp_file::create("gaf-sorter");
        integrated_test(n, 10, 0.05, params);
        temp_file::remove(params.gbwt_file);
    }
}

//------------------------------------------------------------------------------

} // namespace unittest

} // namespace vg

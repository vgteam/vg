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

    GAFInfo(size_t id, size_t nodes, bool with_gbwt_offset) :
        id(id),
        min_node(std::numeric_limits<std::uint32_t>::max()), max_node(0),
        first_node(std::numeric_limits<std::uint32_t>::max()), gbwt_offset(std::numeric_limits<std::uint32_t>::max()) {

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
                return (static_cast<std::uint64_t>(this->first_node) << 32) | this->gbwt_offset;
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

bool is_sorted(const std::string& filename, bool raw_gaf, size_t lines, GAFSorterRecord::key_type key_type) {
    std::ifstream in(filename, std::ios::binary);

    size_t line_num = 0;
    GAFSorterRecord previous;
    while (line_num < lines) {
        GAFSorterRecord record;
        if (raw_gaf) {
            if (!record.read_line(in, key_type)) {
                return false;
            }
        } else {
            if (!record.deserialize(in)) {
                return false;
            }
        }
        if (line_num > 0 && previous.key > record.key) {
            return false;
        }
        previous = record;
        line_num++;
    }

    // There should not be any additional data.
    char c;
    if (in.get(c)) {
        return false;
    }

    return true;
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
        auto lines = generate_gaf(n, 10, 0.05);
        std::string filename = temp_file::create("gaf-sorter");
        GAFSorterFile output(filename);
        sort_gaf_lines(std::move(lines), GAFSorterRecord::key_node_interval, output);
        REQUIRE(is_sorted(filename, true, n, GAFSorterRecord::key_node_interval));
        temp_file::remove(filename);
    }

    SECTION("empty GAF output") {
        size_t n = 0;
        auto lines = generate_gaf(n, 10, 0.05);
        std::string filename = temp_file::create("gaf-sorter");
        GAFSorterFile output(filename);
        sort_gaf_lines(std::move(lines), GAFSorterRecord::key_node_interval, output);
        REQUIRE(is_sorted(filename, true, n, GAFSorterRecord::key_node_interval));
        temp_file::remove(filename);
    }

    SECTION("record output") {
        size_t n = 1234;
        auto lines = generate_gaf(n, 10, 0.05);
        GAFSorterFile output;
        sort_gaf_lines(std::move(lines), GAFSorterRecord::key_node_interval, output);
        REQUIRE(is_sorted(output.name, false, n, GAFSorterRecord::key_node_interval));
        output.remove_temporary();
    }

    SECTION("empty record output") {
        size_t n = 0;
        auto lines = generate_gaf(n, 10, 0.05);
        GAFSorterFile output;
        sort_gaf_lines(std::move(lines), GAFSorterRecord::key_node_interval, output);
        REQUIRE(is_sorted(output.name, false, n, GAFSorterRecord::key_node_interval));
        output.remove_temporary();
    }
}

//------------------------------------------------------------------------------

// FIXME merge, integration

//------------------------------------------------------------------------------

} // namespace unittest

} // namespace vg

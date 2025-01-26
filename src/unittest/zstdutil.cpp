/** \file
 *
 * Unit tests for zstdutil.cpp, which implements wrappers for Zstandard compression.
 */

#include "../zstdutil.hpp"
#include "../utility.hpp"

#include <random>
#include <sstream>

#include "catch.hpp"

namespace vg {

namespace unittest {

//------------------------------------------------------------------------------

namespace {

std::vector<std::string> words {
    "GATTACA", "CAT", "GATTA", "GAGA", "TAG", "GATT"
};

std::string generate_data(size_t num_words, size_t seed) {
    std::string data;
    std::mt19937 rng(num_words ^ seed);
    for (size_t i = 0; i < num_words; i++) {
        size_t index = rng() % words.size();
        data.insert(data.end(), words[index].begin(), words[index].end());
    }
    return data;
}

std::string compress_to_string(const std::string& data) {
    std::stringstream compressed_stream;
    zstd_compress_buf buffer(compressed_stream.rdbuf());
    buffer.sputn(const_cast<char*>(data.data()), data.size());
    buffer.pubsync();
    return compressed_stream.str();
}

void compress_to_file(const std::string& data, const std::string& filename) {
    zstd_ofstream out(filename);
    out.write(data.data(), data.size());
}

std::string decompress_string(const std::string& compressed, size_t expected_size) {
    std::string decompressed;
    decompressed.resize(expected_size);
    std::stringstream compressed_stream(compressed);
    zstd_decompress_buf buffer(compressed_stream.rdbuf());
    size_t n = buffer.sgetn(const_cast<char*>(decompressed.data()), decompressed.size());
    REQUIRE(n == expected_size);
    REQUIRE(buffer.sgetc() == std::char_traits<char>::eof());
    return decompressed;
}

std::string decompress_file(const std::string& filename, size_t expected_size) {
    std::string decompressed;
    decompressed.resize(expected_size);
    zstd_ifstream in(filename);
    in.read(&decompressed[0], decompressed.size());
    REQUIRE(in);
    REQUIRE(in.peek() == std::char_traits<char>::eof());
    return decompressed;
}

} // anonymous namespace

//------------------------------------------------------------------------------

TEST_CASE("Compression with stream buffers", "[zstdutil]") {
    SECTION("empty string") {
        std::string compressed = compress_to_string("");
        std::string decompressed = decompress_string(compressed, 0);
        REQUIRE(decompressed.empty());
    }

    SECTION("random words") {
        for (size_t i = 0; i < 10; i++) {
            std::string data = generate_data(1000, i);
            std::string compressed = compress_to_string(data);
            std::string decompressed = decompress_string(compressed, data.size());
            REQUIRE(decompressed == data);
        }
    }

    SECTION("large instance") {
        for (size_t i = 0; i < 3; i++) {
            std::string data = generate_data(1000000, i);
            std::string compressed = compress_to_string(data);
            std::string decompressed = decompress_string(compressed, data.size());
            REQUIRE(decompressed == data);
        }
    }
}

TEST_CASE("Compression to files", "[zstdutil]") {
    SECTION("empty string") {
        std::string filename = temp_file::create("zstdutil");
        compress_to_file("", filename);
        std::string decompressed = decompress_file(filename, 0);
        REQUIRE(decompressed.empty());
        temp_file::remove(filename);
    }

    SECTION("random words") {
        for (size_t i = 0; i < 10; i++) {
            std::string data = generate_data(1020, i);
            std::string filename = temp_file::create("zstdutil");
            compress_to_file(data, filename);
            std::string decompressed = decompress_file(filename, data.size());
            REQUIRE(decompressed == data);
            temp_file::remove(filename);
        }
    }

    SECTION("large instance") {
        for (size_t i = 0; i < 3; i++) {
            std::string data = generate_data(1000020, i);
            std::string filename = temp_file::create("zstdutil");
            compress_to_file(data, filename);
            std::string decompressed = decompress_file(filename, data.size());
            REQUIRE(decompressed == data);
            temp_file::remove(filename);
        }
    }
}

//------------------------------------------------------------------------------

} // namespace unittest

} // namespace vg

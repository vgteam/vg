#ifndef VG_GFA_NUMERIC_IDS_HPP
#define VG_GFA_NUMERIC_IDS_HPP

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <string_view>
#include <vector>

namespace vg { namespace algorithms {

// Optional index for canonical positive decimal GFA names. Holes remain absent.
// The bitmap is capped at 64 MiB; other names/ranges use the original name map.
class NumericGfaIds {
    bool enabled_;
    std::vector<uint64_t> words;
    size_t count = 0;
public:
    static constexpr uint64_t max_id = (uint64_t(1) << 29) - 1;
    explicit NumericGfaIds(bool enabled = false) : enabled_(enabled) {}
    bool active() const { return enabled_; }
    bool can_store(uint64_t id) const {
        // A sparse outlier should not allocate a large bitmap for a tiny GFA.
        return id <= max_id && id <= std::max(uint64_t(1) << 20, uint64_t(64) * (count + 1));
    }
    static int64_t parse(std::string_view name) {
        if (name.empty() || name.front() < '1' || name.front() > '9') return 0;
        int64_t id = 0;
        for (char c : name) {
            if (c < '0' || c > '9' ||
                id > (std::numeric_limits<int64_t>::max() - (c - '0')) / 10) return 0;
            id = id * 10 + (c - '0');
        }
        return id;
    }
    bool contains(uint64_t id) const {
        return id && id / 64 < words.size() && (words[id / 64] & (uint64_t(1) << (id % 64)));
    }
    bool insert(uint64_t id) {
        size_t required = id / 64 + 1;
        if (required > words.capacity()) {
            size_t capacity = words.capacity() ? words.capacity() : 1;
            while (capacity < required) capacity *= 2;
            words.reserve(capacity);
        }
        if (required > words.size()) words.resize(required, 0);
        uint64_t& word = words[id / 64];
        uint64_t mask = uint64_t(1) << (id % 64);
        bool fresh = !(word & mask);
        word |= mask;
        count += fresh;
        return fresh;
    }
    template<class NameMap>
    void materialize(NameMap& names) {
        for (size_t i = 0; i < words.size(); ++i) {
            uint64_t word = words[i];
            while (word) {
                unsigned bit = __builtin_ctzll(word);
                int64_t id = int64_t(64 * i + bit);
                names.emplace(std::to_string(id), id);
                word &= word - 1;
            }
        }
        std::vector<uint64_t>().swap(words);
        enabled_ = false;
    }
};

}}
#endif

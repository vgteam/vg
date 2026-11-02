#ifndef VG_REFERENCE_PATH_TABLE_HPP_INCLUDED
#define VG_REFERENCE_PATH_TABLE_HPP_INCLUDED

#include <cstddef>
#include <string>
#include <unordered_map>
#include <vector>

namespace vg::multipass {

using namespace std;

/**
 * A directory of the paths vg call places records on, with each one's offset and ploidy.
 *
 * A path's offset is added to every position along it, to place the path on its contig. Its
 * ploidy is the number of copies of the contig the sample has, where no region BED says
 * otherwise.
 *
 * The table does not change once built, so threads can share it.
 */
class ReferencePathTable {
public:
    /// Make an empty table.
    ReferencePathTable() = default;

    /**
     * Make the table of the reference paths `names`.
     *
     * Where a name appears more than once, its last appearance gives its offset and ploidy.
     *
     * @param offsets The offset of each path in `names`, in the same order. A path past the end
     * gets 0.
     * @param ploidies The ploidy of each path in `names`, in the same order. A path past the end
     * gets 2.
     */
    ReferencePathTable(const vector<string>& names, const vector<size_t>& offsets,
                       const vector<int>& ploidies) {
        for (size_t i = 0; i < names.size(); ++i) {
            paths[names[i]] = Path{.offset = i < offsets.size() ? offsets[i] : 0,
                                   .ploidy = i < ploidies.size() ? ploidies[i] : 2};
        }
    }

    /// Return true if the table has no paths.
    bool empty() const { return paths.empty(); }

    /// Return true if `name` is a reference path.
    bool contains(const string& name) const { return paths.contains(name); }

    /// Get the offset of the reference path `name`, or 0 if it is not one.
    size_t offset_of(const string& name) const {
        const auto found = paths.find(name);
        return found != paths.end() ? found->second.offset : 0;
    }

    /// Get the ploidy of the reference path `name`, or 0 if it is not one.
    int ploidy_of(const string& name) const {
        const auto found = paths.find(name);
        return found != paths.end() ? found->second.ploidy : 0;
    }

private:
    /// What the table knows of one reference path.
    struct Path {
        size_t offset = 0;
        int ploidy = 0;
    };
    /// The reference paths, by name.
    unordered_map<string, Path> paths;
};

}

#endif

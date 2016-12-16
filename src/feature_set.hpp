#ifndef VG_FEATURE_SET_H
#define VG_FEATURE_SET_H

/** \file
 *
 * Defines a FeatureSet class that can read and write BED files, and that can
 * keep BED features up to date as paths are edited.
 */

#include <string>
#include <vector>
#include <map>
#include <iostream>
 
namespace vg {

using namespace std;

/**
 * Stores a bunch of features defined on paths, as would be found in a BED file,
 * and listens for messages describing edits to the paths that the features are
 * on. Edit operations are of the form "On path X, range Y to Z has been
 * replaced with a segment of length W". Updates the positions and boundaries of
 * features appropriately.
 *
 * Originally designed to allow BED files to be carried through "vg simplify",
 * but could be used for other annotation liftover tasks.
 */
class FeatureSet {

public:
    /**
     * Represents a Feature that occurs on a path between two inclusive
     * coordinates. Carries along all the extra feature data that BED files
     * store.
     */
    struct Feature {
        /// What Path is the feature on?
        string path_name;
        /// What is the first base's position on that path, inclusive?
        size_t first;
        /// What is the last base's position on the path, inclusive?
        size_t last;
        
        /// What extra BED data should we bring along with this feature?
        /// TODO: BED blocks are not parsed and updated, but they should be.
        vector<string> extra_data;
    };
    
    
    /**
     * Read features from the given BED stream. Adds them to the collection of
     * loaded features.
     */
    void load_bed(istream& in);
    
    /**
     * Save features to the given BED stream.
     */
    void save_bed(ostream& out) const;
    
    /**
     * Notify the FeatureSet that, at the given path, from the given start
     * position, the given number of bases bave been replaced with the other
     * given number of bases. Can handle pure inserts, pure deletions, length-
     * preserving substitutions, and general length-changing substitutions.
     *
     * Updates the contained features that need to change.
     *
     * Note that this is currently O(n) in features on the path. TODO: come up
     * with a truly efficient algorithm to back this using some kind of skip
     * list or rope or something.
     */
    void on_path_edit(string path, size_t start, size_t old_length, size_t new_length);

private:
    /// Stores all the loaded features by path name
    map<string, vector<Feature>> features;

};

}
 
#endif

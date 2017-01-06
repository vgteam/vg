#include "feature_set.hpp"

#include <sstream>

namespace vg {

using namespace std;

void FeatureSet::load_bed(istream& in) {
    // We want to read the BED line by line
    string line;
    while (getline(in, line)) {
    
        // Turn each line into its own stream
        stringstream line_stream(line);
        
        // This line will become a feature
        Feature feature;
        
        // Read the path name
        getline(line_stream, feature.path_name, '\t');
        
        // And the bounds
        line_stream >> feature.first;
        line_stream >> feature.last;
        
        // And the feature name
        line_stream >> feature.feature_name;
        
#ifdef debug
        cerr << "Feature " << feature.feature_name << " from " << feature.first
            << " through " << feature.last << " on " << feature.path_name << endl;
#endif
        
        // TODO: extra data
        
        features[feature.path_name].push_back(feature);
        
    }
}

void FeatureSet::save_bed(ostream& out) const {
    for (auto& kv : features) {
        // For all the contigs
        for (auto& feature : kv.second) {
            // For all the features, dump each one
            out << feature.path_name << "\t" << feature.first << "\t" << feature.last << "\t" << feature.feature_name << endl;
        }
    }
}

void FeatureSet::on_path_edit(const string& path, size_t start, size_t old_length, size_t new_length) {
#ifdef debug
    cerr << "Edit at " << path << " " << start << " from length " << old_length << " to length " << new_length << endl;
#endif

    // Grab the path vector so we can erase if necessary
    auto& path_features = features[path];
    bool erase = false;
    for (auto i = path_features.begin(); i != path_features.end(); i = (erase ? path_features.erase(i) : ++i)) {
        // Go through all the features on the path
        auto& feature = *i;
        
        // Set to true to erase the feature we are on
        erase = false;
        
        if (feature.last < start) {
            // If the feature ends before the start of the edit, do nothing
            continue;
        }
        
        if (feature.first <= start) {
            // Else if it starts at or before the start of the edit
        
            
        
            if (feature.last + 1 < start + old_length) {
                // If it ends before the end of the edit
                
                if (feature.first < start) {
                    // And actually started before the edit
                    
                    // Do some weird interpolation on its end
                    // TODO: interpolate end
                    // Clip for now
                    feature.last = start - 1;
                
#ifdef debug
                    cerr << "\tRight clip feature " << feature.feature_name << " to "
                        << feature.first << " - " << feature.last << endl;
#endif
                } else {
                    // It ends before the edit ends but started at the start of
                    // the edit. We ought to still interpolate, but we're just
                    // going to delete it. TODO: interpolate end.
                    
#ifdef debug
                    cerr << "\tDelete feature " << feature.feature_name << endl;
#endif
                    erase = true;
                    continue;
                }
            } else {
                // If it ends at or after the end of the edit, shift its end up or down by the length difference
                feature.last += ((int64_t) new_length - (int64_t) old_length);
                
#ifdef debug
                cerr << "\tShift end of feature " << feature.feature_name << " by "
                    << ((int64_t) new_length - (int64_t) old_length) << endl;
#endif
            }
            
        } else if (feature.first < start + old_length) {
            // Else if it starts after the start of the edit and before or at the end of the edit
            
            if (feature.last + 1 >= start + old_length) {
                // If it ends after the end of the edit, do some weird
                // interpolation on the start, and shift the end up or down by
                // the length difference. TODO: interpolate start
                
                // Clip for now
                feature.first = start + new_length;
                
                // Adjust end
                feature.last += ((int64_t) new_length - (int64_t) old_length);
                
#ifdef debug
                cerr << "\tLeft clip and shift feature " << feature.feature_name << " to "
                    << feature.first << " - " << feature.last << endl;
#endif
            } else {
                // If it ends at or before the end of the edit, do weird interpolation on the start and end (or just delete it)
                // TODO: interpolate start and end
                // Clip for now
#ifdef debug
                cerr << "\tDelete feature " << feature.feature_name << endl;
#endif
                erase = true;
                continue;
            }
        } else {
            feature.first += ((int64_t) new_length - (int64_t) old_length);
            feature.last += ((int64_t) new_length - (int64_t) old_length);
#ifdef debug
            // Else if it starts after the end of the edit, shift its start and end up or down by the length difference
            cerr << "\tShift feature " << feature.feature_name << " by " << ((int64_t) new_length - (int64_t) old_length) << endl;
#endif
        }
    }
}

const vector<FeatureSet::Feature>& FeatureSet::get_features(const string& path) const {
    return features.at(path);
}

}


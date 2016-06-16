#ifndef VG_TRANSLATOR_H
#define VG_TRANSLATOR_H
// translator.hpp: defines the Translator, which maps paths on an augmented graph
// into the base graph defined by a set of Translations from the augmented graph

#include <iostream>
#include <algorithm>
#include <functional>
#include <set>
#include <vector>
#include <list>
#include "vg.pb.h"
#include "vg.hpp"
#include "hash_map.hpp"
#include "utility.hpp"
#include "types.hpp"

namespace vg {

using namespace std;

/**
 * Class to map paths into a base graph found via a set of Translations
 */
class Translator {
public:

    vector<Translation> translations;
    map<pos_t, Translation*> pos_to_trans;
    Translator(const vector<Translation>& trans);
    void load_translations(const vector<Translation>& trans);
    Position translate(const Position& position);
    Path translate(const Path& path);
    
};



}


#endif

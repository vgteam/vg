#ifndef VG_TRANSLATOR_HPP_INCLUDED
#define VG_TRANSLATOR_HPP_INCLUDED
/// \file translator.hpp
/// Defines the Translator, which maps paths on an augmented graph into the base
/// graph defined by a set of Translations from the augmented graph

#include <iostream>
#include <algorithm>
#include <functional>
#include <set>
#include <vector>
#include <list>
#include <sstream>
#include <exception>
#include <vg/vg.pb.h>
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
    Translator(void);
    Translator(istream& in);
    Translator(const vector<Translation>& trans);
    void load(const vector<Translation>& trans);
    void build_position_table(void);
    Translation get_translation(const Position& position);
    bool has_translation(const Position& position, bool ignore_strand = true);
    Position translate(const Position& position);
    Position translate(const Position& position, const Translation& translation);
    Edge translate(const Edge& edge);
    Mapping translate(const Mapping& mapping);
    Path translate(const Path& path);
    Alignment translate(const Alignment& aln);
    Locus translate(const Locus& locus);
    Translation overlay(const Translation& trans);
};

bool is_match(const Translation& translation);

}


#endif

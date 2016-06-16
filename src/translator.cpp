#include "translator.hpp"

namespace vg {

using namespace std;

Translator::Translator(const vector<Translation>& trans) {
    load_translations(trans);
}

void Translator::load_translations(const vector<Translation>& trans) {
    translations = trans;
    //map<pos_t, Translation*> pos_to_trans;
    for (auto& t : translations) {
        // map from the new positions to the corresponding translations
        pos_to_trans[make_pos_t(t.to().mapping(0).position())] = &t;
    }
}

Position Translator::translate(const Position& position) {
    auto pos = make_pos_t(position);
    auto t = pos_to_trans.lower_bound(pos);
    return t->second->from().mapping(0).position();
}

Path Translator::translate(const Path& path) {
    
}

}

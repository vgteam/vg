#include "translator.hpp"
#include "stream.hpp"

namespace vg {

using namespace std;

Translator::Translator(istream& in) {
    function<void(Translation&)> lambda = [&](Translation& trans) {
        translations.push_back(trans);
    };
    stream::for_each(in, lambda);
    build_position_table();
}

Translator::Translator(const vector<Translation>& trans) {
    translations = trans;
    build_position_table();
}

void Translator::build_position_table(void) {
    //map<pos_t, Translation*> pos_to_trans;
    for (auto& t : translations) {
        // map from the new positions to the corresponding translations
        pos_to_trans[make_pos_t(t.to().mapping(0).position())] = &t;
    }
}

Translation Translator::get_translation(const Position& position) {
    auto pos = make_pos_t(position);
    // check that the node is in the translation
    pos_t node_only = pos;
    get_offset(node_only) = 0;
    get_is_rev(node_only) = false;
    Translation translation;
    if (pos_to_trans.find(node_only) == pos_to_trans.end()) {
        cerr << "WARNING: node " << id(pos) << " is not in the translation table" << endl;
    } else {
        auto t = pos_to_trans.upper_bound(pos); --t;
        translation = *t->second;
    }
    return translation;
}

Position Translator::translate(const Position& position) {
    return translate(position, get_translation(position));
}

Position Translator::translate(const Position& position, const Translation& translation) {
    // what kind of translation is it?
    if (is_match(translation)) {
        if (position.offset() >= mapping_from_length(translation.to().mapping(0))) {
            cerr << "ERROR: to-position offset is greater than translation length "
                 << mapping_from_length(translation.to().mapping(0)) << endl;
            exit(1);
        }
        // exact match; local coordinate space is identical
        Position from_pos = translation.from().mapping(0).position();
        from_pos.set_offset(from_pos.offset()+position.offset());
        return from_pos;
    } else {
        // novel sequence
        return translation.from().mapping(0).position();
    }
}

// TODO MAKE THERE BE A REVERSE TRANSLATIONS OBJECT
// AND STORE THEM REVERSE TRANSLATIONS IN IT

Mapping Translator::translate(const Mapping& mapping) {
    Mapping translated = mapping;
    if (!mapping.has_position()) return mapping;
    Translation translation = get_translation(mapping.position());
    *translated.mutable_position() = translate(mapping.position(), translation);
    if (is_match(translation)) {
        return translated;
    } else {
        // subset the translated path mapping.from_length() bases in or out
        // depending on its offset
        //const string mapping_sequence(const Mapping& m, const string& node_seq);
        if (mapping_is_match(mapping)) {
            if (mapping.position().offset() == 0) {
                // we can take the position from the translation
                // so subset from the beginning of the translation
                translated =
                    cut_mapping(
                        cut_mapping(translation.from().mapping(0),
                                    mapping_from_length(mapping)).first,
                        mapping_to_length(mapping)).first;
            } else if (mapping.rank() == 1) {
                // if we're the first but we don't have a position, it's possible that we are going to reach
                // the old reference graph in the next step
                translated =
                    cut_mapping(
                        cut_mapping(translation.from().mapping(0),
                                    mapping_from_length(mapping)).second,
                        mapping_to_length(mapping)).second;
                //auto seq = mapping_sequence(mapping, translated.edit(0).sequence());
                //translated.mutable_edit(0)->set_sequence(mapping_sequence(mapping, translated.edit(0).sequence()));
                translated.clear_position();
            } else {
                // clear the position as we'd be unmapped in the from graph
                translated.clear_position();
            }
        } else {
            if (mapping.position().offset() == 0) {
                // we can take the position from the translation
                // so subset from the beginning of the translation
                translated =
                    cut_mapping(translation.from().mapping(0),
                                mapping_from_length(mapping)).first;
                auto seq = translation.from().mapping(0).edit(0).sequence();
                translated.clear_edit();
                auto edit = translated.add_edit();
                edit->set_sequence(mapping_sequence(mapping, seq));
                edit->set_to_length(edit->sequence().size());
            } else if (mapping.rank() == 1) {
                // if we're the first but we don't have a position, it's possible that we are going to reach
                // the old reference graph in the next step
                translated =
                    cut_mapping(translation.from().mapping(0),
                                mapping_from_length(mapping)).second;
                auto seq = translation.from().mapping(0).edit(0).sequence();
                translated.clear_edit();
                auto edit = translated.add_edit();
                edit->set_sequence(mapping_sequence(mapping, seq));
                edit->set_to_length(edit->sequence().size());
                translated.clear_position();
            } else {
                // clear the position as we'd be unmapped in the from graph
                translated.clear_position();
            }
        }
        return translated;
    }
}

Path Translator::translate(const Path& path) {
    Path result;
    for (int i = 0; i < path.mapping_size(); ++i) {
        *result.add_mapping() = translate(path.mapping(i));
    }
    return result;
}

Alignment Translator::translate(const Alignment& aln) {
    Alignment result = aln;
    *result.mutable_path() = translate(aln.path());
    return result;
}

bool is_match(const Translation& translation) {
    return
        path_is_simple_match(translation.from())
        && path_is_simple_match(translation.to())
        && path_from_length(translation.from()) == path_to_length(translation.from())
        && path_from_length(translation.to()) == path_to_length(translation.to())
        && path_to_length(translation.from()) == path_to_length(translation.to());
}

}

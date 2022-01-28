/**
 * \file back_translate.cpp
 */

#include "back_translate.hpp"
#include "../path.hpp"

namespace vg {
namespace algorithms {

using namespace std;

/**
 * Erase the items at the given indices in the given Protobuf repeated pointer field.
 */
template<typename RepeatedPtrField>
static void erase_at(RepeatedPtrField* field, const vector<size_t>& sorted_indices_to_remove) {
    if (sorted_indices_to_remove.empty()) {
        // We don't need to do anything.
        return;
    }
    
    // This is how many slots ahead we are reading to write the current slot.
    // It's also (1+) the cursor in sorted_indices_to_remove
    size_t slots_ahead = 1;
    for (size_t i = sorted_indices_to_remove[slots_ahead - 1]; i + sorted_indices_to_remove.size() < field->size(); i++) {
        // We're at a slot that will need to be filled in the output.
        while (slots_ahead < sorted_indices_to_remove.size() && i + slots_ahead == sorted_indices_to_remove[slots_ahead]) {
            // Slide copy source ahead until it stops being the next thing to skip.
            slots_ahead++;
        }
        // Overwrite the item here with the item that is supposed to be here.
        // We use a swap so we don't have to actually copy any item data, just pointers.
        field->SwapElements(i, i + slots_ahead);
    }
    // Now we've bumped all the unused items to the end so delete them.
    field->DeleteSubrange(field->size() - sorted_indices_to_remove.size(), sorted_indices_to_remove.size());
}

void back_translate_in_place(const NamedNodeBackTranslation* translation, Path& path) {
    // For combining mappings that can combine, we need to track where the previous mapping ended.
    nid_t prev_segment_number = numeric_limits<nid_t>::max();
    bool prev_segment_is_reverse = false;
    size_t prev_segment_offset = numeric_limits<size_t>::max();
    // And also what mapping it is or has been glommed into
    Mapping* prev_mapping = nullptr;
    
    // When we glom mappings into other mappings we steal their edits and put
    // their indices in this list.
    vector<size_t> mapping_indices_to_remove;
    
    for (size_t i = 0; i < path.mapping_size(); i++) {
        // For each Mapping
        Mapping* mapping = path.mutable_mapping(i);
        
        // Determine its range
        oriented_node_range_t source_range(mapping->position().node_id(), mapping->position().is_reverse(),
                                           mapping->position().offset(), from_length(*mapping)); 
        
        // Translate it
        vector<oriented_node_range_t> translated = translation->translate_back(source_range);
        
        if (translated.size() != 1) {
            // TODO: Implement translations that split graph nodes into multiple segments.
            throw std::runtime_error("Translated range on node " + to_string(get<0>(source_range)) +
                                     " to " + to_string(translated.size()) + 
                                     " named segment ranges, but complex translations like this are not yet implemented");
        }
        
        auto& translated_range = translated[0];
        if (get<1>(translated_range) != get<1>(source_range)) {
            // TODO: Implement translations that flip orientations.
            throw std::runtime_error("Translated range on node " + to_string(get<0>(source_range)) +
                                     " ended up on the opposite strand; complex translations like this are not yet implemented");
        }
        
        // Change over to the named sequence and the offset there.
        mapping->mutable_position()->clear_node_id();
        mapping->mutable_position()->set_name(translation->get_back_graph_node_name(get<0>(translated_range)));
        mapping->mutable_position()->set_offset(get<2>(translated_range));
        
        if (i == 0 || get<0>(translated_range) != prev_segment_number || get<1>(translated_range) != prev_segment_is_reverse || get<2>(translated_range) != prev_segment_offset) {
            // We have done a transition that isn't just abutting in a
            // segment. We assume anything we want to represent as a
            // deletion is already represented as a deletion, so we
            // preserve all jumps as jumps.
            
            // Just move to this part of this segment, and then advance by
            // the length of the piece we translated.
            prev_segment_number = get<0>(translated_range);
            prev_segment_is_reverse = get<1>(translated_range);
            prev_segment_offset = get<2>(translated_range) + get<3>(translated_range);
            
            // We are now the prev mapping to glom into
            prev_mapping = mapping;
        } else {
            // We abut the previous Mapping. So we should be able to glom
            // all our edits into it, and then stop existing.
            for (size_t j = 0; j < mapping->edit_size(); j++) {
                // Glom each edit into the previous mapping
                
                if (j == 0 && prev_mapping->edit_size() > 0 && edits_are_compatible(prev_mapping->edit(prev_mapping->edit_size() - 1), mapping->edit(j))) {
                    // We assume our edits are all merged up if they can be. So
                    // we only have to consider merging the first of our edits
                    // into the last edit of the previous mapping. Turns out
                    // they can merge.
                    merge_edits_in_place(*prev_mapping->mutable_edit(prev_mapping->edit_size() - 1), mapping->edit(j));
                } else {
                    // This edit has to become a new edit in the previous mapping.
                    *prev_mapping->add_edit() = std::move(*mapping->mutable_edit(j));
                }
            }
            
            // Now we need to get rid of this mapping.
            // Leave prev_mapping pointing where it is.
            // And remember we don't want this mapping anymore.
            mapping_indices_to_remove.push_back(i);
            
            // Advance along the segment
            prev_segment_offset += get<3>(translated_range);
        }
    }
    
    // We need to batch-erase all these indices.
    erase_at(path.mutable_mapping(), mapping_indices_to_remove);
}



void back_translate_in_place(const NamedNodeBackTranslation* translation, Snarl& snarl) {
    // To translate a snarl, you translate its bounding visits
    back_translate_in_place(translation, *snarl.mutable_start());
    back_translate_in_place(translation, *snarl.mutable_end());
}

void back_translate_in_place(const NamedNodeBackTranslation* translation, SnarlTraversal& traversal) {
    vector<size_t> visit_indices_to_remove;
    for (size_t i = 0; i < traversal.visit_size(); i++) {
        // Translate every visit
        Visit& here = *(traversal.mutable_visit(i));
        back_translate_in_place(translation, here);
        if (i > 0) {
            const Visit& prev = traversal.visit(i - 1);
            if (!here.has_snarl() && !prev.has_snarl() && here.name() == prev.name() && here.backward() == prev.backward()) {
                // These visits can coalesce because they are to the same segment.
                visit_indices_to_remove.push_back(i);
            }
        }
    }
    
    // Get rid of visits that coalesced away
    erase_at(traversal.mutable_visit(), visit_indices_to_remove);
}

void back_translate_in_place(const NamedNodeBackTranslation* translation, Visit& visit) {
    if (visit.has_snarl()) {
        // A visit with a snarl is translated by translating its snarl
        back_translate_in_place(translation, *visit.mutable_snarl());
    } else {
        // Otherwise, translate its node ID to a segment name. Use made-up boundaries on the node.
        // TODO: Can we have an easy way to say "whole node"?
        oriented_node_range_t source_range(visit.node_id(), visit.backward(), 0, 1);
        
        // Translate it
        vector<oriented_node_range_t> translated = translation->translate_back(source_range);
        
        if (translated.size() != 1) {
            // TODO: Implement translations that split graph nodes into multiple segments.
            throw std::runtime_error("Translated range on node " + to_string(get<0>(source_range)) +
                                     " to " + to_string(translated.size()) + 
                                     " named segment ranges, but complex translations like this are not yet implemented");
        }
        
        auto& translated_range = translated[0];
        if (get<1>(translated_range) != get<1>(source_range)) {
            // TODO: Implement translations that flip orientations.
            throw std::runtime_error("Translated range on node " + to_string(get<0>(source_range)) +
                                     " ended up on the opposite strand; complex translations like this are not yet implemented");
        }
        
        // Save the change to the visit.
        visit.clear_node_id();
        visit.set_name(translation->get_back_graph_node_name(get<0>(translated_range)));
        // Ignore the offset and length info.
    }
}

}
}

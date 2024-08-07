/// \file annotation.hpp
/// utilities for working with Struct-formatted annotations on Protobuf objects


#ifndef VG_ANNOTATION_HPP_INCLUDED
#define VG_ANNOTATION_HPP_INCLUDED

#include <vector>
#include <string>
#include <type_traits>
#include <functional>
#include <limits>
#include <cmath>
#include <sstream>

#include <google/protobuf/struct.pb.h>

namespace vg {

using namespace std;

////////////////////////////////////////////////////////////////////////
// API
////////////////////////////////////////////////////////////////////////

/// Returns true if the Protobuf object has an annotation with this name
template<typename Annotated>
bool has_annotation(const Annotated& annotated, const string& name);

/// Get the annotation with the given name and return it.
/// If not present, returns the Protobuf default value for the annotation type.
/// The value may be a primitive type or an entire Protobuf object.
/// It is undefined behavior to read a value out into a different type than it was stored with.
template<typename AnnotationType, typename Annotated>
AnnotationType get_annotation(const Annotated& annotated, const string& name);

/// Get the annotation with the given name and return it.
/// If not present, returns the Protobuf default value for the annotation type.
/// The value may be a primitive type or an entire Protobuf object.
/// It is undefined behavior to read a value out into a different type than it was stored with.
template<typename AnnotationType, typename Annotated>
AnnotationType get_annotation(Annotated* annotated, const string& name);

// We implement the mutators vor pointers (to match the Protobuf
// mutable_whatever()) as well as for references (because I kept forgetting to
// use the pointers)

/// Set the annotation with the given name to the given value. The value may be
/// a primitive type or a vector of a primitive type.
template<typename AnnotationType, typename Annotated>
void set_annotation(Annotated* annotated, const string& name, const AnnotationType& annotation);

/// Set the annotation with the given name to the given value.
/// The value may be a primitive type or a vector of a primitive type.
template<typename AnnotationType, typename Annotated>
void set_annotation(Annotated& annotated, const string& name, const AnnotationType& annotation);

/// Set a pair of annotations to compactly express the values in the given
/// vector which contains many repeated values. The values will be sorted in place.
template<typename AnnotatonValueType, typename Annotated>
void set_compressed_annotation(Annotated* annotated, const string& base_name, std::vector<AnnotatonValueType> annotation);

/// Set a pair of annotations to compactly express the values in the given
/// vector which contains many repeated values. The values will be sorted in place.
template<typename AnnotatonValueType, typename Annotated>
void set_compressed_annotation(Annotated& annotated, const string& base_name, std::vector<AnnotatonValueType> annotation);

/// Clear the annotation with the given name.
template<typename Annotated>
void clear_annotation(Annotated* annotated, const string& name);

/// Clear the annotation with the given name
template<typename Annotated>
void clear_annotation(Annotated& annotated, const string& name);

/// Apply a lambda to all annotations, except for Struct and ListValue annotations (which cannot
/// be easily typed without exposing ugly Protobuf internals
template<typename Annotated>
void for_each_basic_annotation(const Annotated& annotated,
                               const function<void(const string&)> null_lambda,
                               const function<void(const string&,double)> double_lambda,
                               const function<void(const string&,bool)> bool_lambda,
                               const function<void(const string&,const string&)> string_lambda);

////////////////////////////////////////////////////////////////////////
// Internal Definitions
////////////////////////////////////////////////////////////////////////

/// We define an adapter for things that are annotated to let us get at the
/// annotation struct. It is only defined for the actual types (Alignment,
/// MultipathAlignment) and not pointers to them. This keeps the API overloads
/// that are supposed to be for references to the types from operating on
/// references to pointers to the types instead.
template<typename T, typename Enabled = typename enable_if<!is_pointer<T>::value>::type>
struct Annotation {
    /// Get the immutable annotations Struct
    static const google::protobuf::Struct& get(const T& t);
    /// Get the mutable annotations struct.
    static google::protobuf::Struct* get_mutable(T* t);
    /// Clear all annotations
    void clear(T* t);
};

/// Cast a Protobuf generic Value to any type.
template<typename T>
inline T value_cast(const google::protobuf::Value& value);

/// Cast any type to a generic Protobuf value.
template<typename T>
inline google::protobuf::Value value_cast(const T& wrap);

////////////////////////////////////////////////////////////////////////
// Implementation
////////////////////////////////////////////////////////////////////////

template<typename T, typename Enabled>
const google::protobuf::Struct& Annotation<T, Enabled>::get(const T& t) {
    return t.annotation();
}

template<typename T, typename Enabled>
google::protobuf::Struct* Annotation<T, Enabled>::get_mutable(T* t) {
    return t->mutable_annotation();
}

template<typename T, typename Enabled>
void Annotation<T, Enabled>::clear(T* t) {
    t->clear_annotation();
}

// We define all these value_cast implementations, in both directions

// For Struct we use a pointer so you can tell if it's not really there by having a nullptr.
template<>
inline const google::protobuf::Struct* value_cast<const google::protobuf::Struct*>(const google::protobuf::Value& value) {
    assert(value.kind_case() == google::protobuf::Value::KindCase::kStructValue);
    return &value.struct_value();
}

// For Value we use a pointer so you can tell if it's not really there by having a nullptr.
template<>
inline const google::protobuf::Value* value_cast<const google::protobuf::Value*>(const google::protobuf::Value& value) {
    return &value;
}

template<>
inline bool value_cast<bool>(const google::protobuf::Value& value) {
    assert(value.kind_case() == google::protobuf::Value::KindCase::kBoolValue);
    return value.bool_value();
}

template<>
inline double value_cast<double>(const google::protobuf::Value& value) {
    if (value.kind_case() == google::protobuf::Value::KindCase::kStringValue) {
        // If someone puts in an infinite or NAN double, Protobuf refuses to
        // stringify those, so we do it ourselves. But now they want the double
        // back so we need to undo that.
        if (value.string_value() == "Infinity") {
            return std::numeric_limits<double>::infinity();
        } else if (value.string_value() == "-Infinity") {
            return -std::numeric_limits<double>::infinity();
        } else if (value.string_value() == "NaN") {
            return nan("");
        } else {
            throw std::runtime_error("Cannot understand " + value.string_value() + " as a double.");
        }
    }
    assert(value.kind_case() == google::protobuf::Value::KindCase::kNumberValue);
    return value.number_value();
}

template<>
inline string value_cast<string>(const google::protobuf::Value& value) {
    assert(value.kind_case() == google::protobuf::Value::KindCase::kStringValue);
    return value.string_value();
}

template<>
inline google::protobuf::Value value_cast<bool>(const bool& wrap) {
    google::protobuf::Value to_return;
    to_return.set_bool_value(wrap);
    return to_return;
}

template<>
inline google::protobuf::Value value_cast<double>(const double& wrap) {
    google::protobuf::Value to_return;
    // We need to represent inf and nan values as something else, since Protobuf now refuses to serialize them as anything.
    // Previously it made them "Infinity", "-Infinity" and "NaN", so we do that too.
    if (isinf(wrap)) {
        to_return.set_string_value(wrap > 0 ? "Infinity" : "-Infinity");
    } else if (isnan(wrap)) {
        to_return.set_string_value("NaN");
    } else {
        to_return.set_number_value(wrap);
    }
    return to_return;
}

template<>
inline google::protobuf::Value value_cast<string>(const string& wrap) {
    google::protobuf::Value to_return;
    to_return.set_string_value(wrap);
    return to_return;
}

// Helpers for dumping integral types to double.
// May lose precision for large numbers.

template<>
inline google::protobuf::Value value_cast<size_t>(const size_t& wrap) {
    return value_cast<double>((double) wrap);
}

template<>
inline google::protobuf::Value value_cast<int>(const int& wrap) {
    return value_cast<double>((double) wrap);
}

// We also have implementations for vectors and other push_back-able containers.

template<typename Container>
inline Container value_cast(const google::protobuf::Value& value) {
    assert(value.kind_case() == google::protobuf::Value::KindCase::kListValue);
    Container items;
    for (auto& nested_value : value.list_value().values()) {
        items.push_back(value_cast<typename Container::value_type>(nested_value));
    }
    return items;
}

template<typename Container>
inline google::protobuf::Value value_cast(const Container& wrap) {
    // Make a new list that the Protobuf Value message will eventually own
    google::protobuf::ListValue* list = new google::protobuf::ListValue();
    for (auto& item : wrap) {
        // Put all the items from the container into it
        *list->add_values() = value_cast(item);
    }

    google::protobuf::Value to_return;
    // Hand it off
    to_return.set_allocated_list_value(list);
    return to_return;
}

template<typename Annotated>
bool has_annotation(const Annotated& annotated, const string& name) {
    // Grab the whole annotation struct
    const google::protobuf::Struct& annotation_struct = Annotation<Annotated>::get(annotated);

    const google::protobuf::Struct* here = &annotation_struct;
    const google::protobuf::Value* leaf = nullptr;
    std::string name_part;
    std::istringstream ss(name);
    while (std::getline(ss, name_part, '.')) {
        if (here == nullptr) {
            // Path extends beyond a leaf value
            return false;
        }
        // Look up each dot-separated segment
        auto found = here->fields().find(name_part);
        if (found == here->fields().end()) {
            // This segment isn't present
            return false;
        }
        const google::protobuf::Value& part_value = found->second;
        if (part_value.kind_case() == google::protobuf::Value::KindCase::kStructValue) {
            // Recurse into the struct
            here = &part_value.struct_value();
        } else {
            // Maybe this is the last segment and we found the actual thing?
            here = nullptr;
            leaf = &part_value;
        }
    }
    // If we get here, we ran out of name
    // Return true if there is any value here, even a struct
    return true;
}

// TODO: more value casts for e.g. ints and embedded messages.

template<typename AnnotationType, typename Annotated>
AnnotationType get_annotation(const Annotated& annotated, const string& name) {
    // Grab the whole annotation struct
    const google::protobuf::Struct& annotation_struct = Annotation<Annotated>::get(annotated);

    const google::protobuf::Struct* here = &annotation_struct;
    const google::protobuf::Value* leaf = nullptr;
    std::string name_part;
    std::istringstream ss(name);
    while (std::getline(ss, name_part, '.')) {
        if (here == nullptr) {
            // Path extends beyond a leaf value
            // Return the Proto default value, by value-initializing.
            return AnnotationType();
        }
        // Look up each dot-separated segment.
        // We don't use find because the find interface can't get us references
        // into the Protobuf storage for giving back Value or Struct pointers.
        if (!here->fields().count(name_part)) {
            // This segment isn't present
            // Return the Proto default value, by value-initializing.
            return AnnotationType();
        }
        const google::protobuf::Value& part_value = here->fields().at(name_part);
        if (part_value.kind_case() == google::protobuf::Value::KindCase::kStructValue) {
            // Recurse into the struct
            here = &part_value.struct_value();
            // We might be fetching the whole struct though
            leaf = &part_value;
        } else {
            // Maybe this is the last segment and we found the actual thing?
            here = nullptr;
            leaf = &part_value;
        }
    }
    
    // Pull out the right type from the leaf Value.
    return value_cast<AnnotationType>(*leaf);
}

template<typename AnnotationType, typename Annotated>
inline AnnotationType get_annotation(Annotated* annotated, const string& name) {
    return get_annotation<AnnotationType>(*annotated, name);
}

template<typename AnnotationType, typename Annotated>
void set_annotation(Annotated* annotated, const string& name, const AnnotationType& annotation) {
    // Get ahold of the struct
    google::protobuf::Struct* annotation_struct = Annotation<Annotated>::get_mutable(annotated);

    google::protobuf::Struct* here = annotation_struct;
    google::protobuf::Value* leaf = nullptr;
    std::string name_part;
    std::istringstream ss(name);
    while (std::getline(ss, name_part, '.')) {
        // Look up each dot-separated segment and put a struct there
        leaf = &(*here->mutable_fields())[name_part];
        here = leaf->mutable_struct_value();
    }

    assert(leaf != nullptr);

    // Actually make the last one not a struct but a real leaf value
    here = nullptr;
    *leaf = value_cast(annotation);
}

template<typename AnnotationType, typename Annotated>
inline void set_annotation(Annotated& annotated, const string& name, const AnnotationType& annotation) {
    set_annotation(&annotated, name, annotation);
}

template<typename AnnotatonValueType, typename Annotated>
void set_compressed_annotation(Annotated* annotated, const string& base_name, std::vector<AnnotatonValueType> annotation) {
    // Sort the values
    std::sort(annotation.begin(), annotation.end());
    
    std::vector<AnnotatonValueType> values;
    std::vector<size_t> counts;
    bool duplicates = false;
    for (auto& v : annotation) {
        // Run lenght compress the values
        if (!values.empty() && v == values.back()) {
            counts.back()++;
            duplicates = true;
        } else {
            values.push_back(v);
            counts.push_back(1);
        }
    }

    // Apply two annotations
    set_annotation(annotated, base_name + ".values", values);
    if (duplicates) {
        // Only include the weights if some are not 1
        set_annotation(annotated, base_name + ".weights", counts);
    }
}

template<typename AnnotatonValueType, typename Annotated>
inline void set_compressed_annotation(Annotated& annotated, const string& base_name, std::vector<AnnotatonValueType> annotation) {
    set_compressed_annotation(&annotated, base_name, annotation);
}

template<typename Annotated>
void clear_annotation(Annotated* annotated, const string& name) {
    // Get ahold of the struct
    google::protobuf::Struct* annotation_struct = Annotation<Annotated>::get_mutable(annotated);
    
    google::protobuf::Struct* parent = nullptr;
    google::protobuf::Struct* here = annotation_struct;
    std::string name_part;
    std::string last_part;
    std::istringstream ss(name);
    while (std::getline(ss, name_part, '.')) {
        if (here == nullptr) {
            // Path extends beyond a leaf value
            return;
        }
        // Look up each dot-separated segment
        auto found = here->mutable_fields()->find(name_part);
        if (found == here->mutable_fields()->end()) {
            // This segment isn't present
            return;
        }
        google::protobuf::Value* part_value = &found->second;
        if (part_value->kind_case() == google::protobuf::Value::KindCase::kStructValue) {
            // Recurse into the struct
            parent = here;
            here = part_value->mutable_struct_value();
        } else {
            // Maybe this is the last segment and we found the actual thing?
            parent = here;
            here = nullptr;
        }
        last_part = std::move(name_part);
    }

    if (parent != nullptr) {
        // Clear out that field
        here = nullptr;
        parent->mutable_fields()->erase(last_part);
    }
}

template<typename Annotated>
inline void clear_annotation(Annotated& annotated, const string& name) {
    clear_annotation(&annotated, name);
}

template<typename Annotated>
void for_each_basic_annotation(const Annotated& annotated,
                               const function<void(const string&)> null_lambda,
                               const function<void(const string&,double)> double_lambda,
                               const function<void(const string&,bool)> bool_lambda,
                               const function<void(const string&,const string&)> string_lambda) {
    
    for (auto it = annotated.annotation().fields().begin(), end = annotated.annotation().fields().end(); it != end; ++it) {
        switch (it->second.kind_case()) {
            case google::protobuf::Value::KindCase::kBoolValue:
                bool_lambda(it->first, it->second.bool_value());
                break;
            case google::protobuf::Value::KindCase::kNumberValue:
                double_lambda(it->first, it->second.number_value());
                break;
            case google::protobuf::Value::KindCase::kStringValue:
                string_lambda(it->first, it->second.string_value());
                break;
            case google::protobuf::Value::KindCase::kNullValue:
                null_lambda(it->first);
                break;
            default:
                // TODO: skip ListValue and Struct, how to include?
                break;
        }
    }
}

}



#endif

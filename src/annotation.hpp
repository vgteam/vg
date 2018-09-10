/// \file annotation.hpp
/// utilities for working with Struct-formatted annotations on Protobuf objects


#ifndef VG_ANNOTATION_HPP_INCLUDED
#define VG_ANNOTATION_HPP_INCLUDED

#include <vector>
#include <string>
#include <type_traits>

#include <google/protobuf/struct.pb.h>

namespace vg {

using namespace std;

////////////////////////////////////////////////////////////////////////
// API
////////////////////////////////////////////////////////////////////////

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

/// Clear the annotation with the given name.
template<typename Annotated>
void clear_annotation(Annotated* annotated, const string& name);

/// Clear the annotation with the given name
template<typename Annotated>
void clear_annotation(Annotated& annotated, const string& name);

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
template <typename T>
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

template<>
inline bool value_cast<bool>(const google::protobuf::Value& value) {
    assert(value.kind_case() == google::protobuf::Value::KindCase::kBoolValue);
    return value.bool_value();
}

template<>
inline double value_cast<double>(const google::protobuf::Value& value) {
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
    to_return.set_number_value(wrap);
    return to_return;
}

template<>
inline google::protobuf::Value value_cast<string>(const string& wrap) {
    google::protobuf::Value to_return;
    to_return.set_string_value(wrap);
    return to_return;
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

// TODO: more value casts for e.g. ints and embedded messages.

template<typename AnnotationType, typename Annotated>
inline AnnotationType get_annotation(const Annotated& annotated, const string& name) {
    // Grab the whole annotation struct
    auto annotation_struct = Annotation<Annotated>::get(annotated);
    
    if (!annotation_struct.fields().count(name)) {
        // Nothing is there.
        // Return the Proto default value, by value-initializing.
        return AnnotationType();
    }
    
    // Get the Protobuf Value for this annotation name
    auto value = annotation_struct.fields().at(name);
    
    // Pull out the right type.
    return value_cast<AnnotationType>(value);
}

template<typename AnnotationType, typename Annotated>
inline AnnotationType get_annotation(Annotated* annotated, const string& name) {
    return get_annotation<AnnotationType>(*annotated, name);
}

template<typename AnnotationType, typename Annotated>
inline void set_annotation(Annotated* annotated, const string& name, const AnnotationType& annotation) {
    // Get ahold of the struct
    auto* annotation_struct = Annotation<Annotated>::get_mutable(annotated);
    
    // Set the key to the wrapped value
    (*annotation_struct->mutable_fields())[name] = value_cast(annotation);
}

template<typename AnnotationType, typename Annotated>
inline void set_annotation(Annotated& annotated, const string& name, const AnnotationType& annotation) {
    set_annotation(&annotated, name, annotation);
}

template<typename Annotated>
inline void clear_annotation(Annotated* annotated, const string& name) {
    // Get ahold of the struct
    auto* annotation_struct = Annotation<Annotated>::get_mutable(annotated);
    // Clear out that field
    annotation_struct->mutable_fields()->erase(name);
}

template<typename Annotated>
inline void clear_annotation(Annotated& annotated, const string& name) {
    clear_annotation(&annotated, name);
}

}



#endif

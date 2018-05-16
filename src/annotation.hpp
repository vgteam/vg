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

/// Set the annotation with the given name to the given value.
/// The value may be a primitive type or an entire Protobuf object.
template<typename AnnotationType, typename Annotated>
void set_annotation(Annotated* annotated, const string& name, const AnnotationType& annotation);

/// Clear the annotation with the given name.
template<typename Annotated>
void clear_annotation(Annotated* annotated, const string& name);

////////////////////////////////////////////////////////////////////////
// Internal Definitions
////////////////////////////////////////////////////////////////////////

/// We define an adapter for things that are annotated to let us get at the annotation struct.
template<typename T>
struct Annotation {
    /// Get the immutable annotations Struct
    static const google::protobuf::Struct& get(const T& t);
    /// Get the mutable annotations struct.
    static google::protobuf::Struct* get_mutable(T* t);
    /// Clear all annotations
    void clear(T* t);
};

/// Cast a Protobuf generic Value to any type.
template <typename T, typename Enabled = void>
inline T value_cast(const google::protobuf::Value& value);

/// Cast any type to a generic Protobuf value.
template<typename T,  typename Enabled = void>
inline google::protobuf::Value value_cast(const T& wrap);

////////////////////////////////////////////////////////////////////////
// Implementation
////////////////////////////////////////////////////////////////////////

template<typename T>
const google::protobuf::Struct& Annotation<T>::get(const T& t) {
    return t.annotation();
}

template<typename T>
google::protobuf::Struct* Annotation<T>::get_mutable(T* t) {
    return t->mutable_annotation();
}

template<typename T>
void Annotation<T>::clear(T* t) {
    t->clear_annotation();
}

template<>
inline bool value_cast<bool, void>(const google::protobuf::Value& value) {
    assert(value.kind_case() == google::protobuf::Value::KindCase::kBoolValue);
    return value.bool_value();
}

template<>
inline double value_cast<double, void>(const google::protobuf::Value& value) {
    assert(value.kind_case() == google::protobuf::Value::KindCase::kNumberValue);
    return value.number_value();
}

template<>
inline string value_cast<string, void>(const google::protobuf::Value& value) {
    assert(value.kind_case() == google::protobuf::Value::KindCase::kStringValue);
    return value.string_value();
}

template<>
inline google::protobuf::Value value_cast<bool, void>(const bool& wrap) {
    google::protobuf::Value to_return;
    to_return.set_bool_value(wrap);
    return to_return;
}

template<>
inline google::protobuf::Value value_cast<double, void>(const double& wrap) {
    google::protobuf::Value to_return;
    to_return.set_number_value(wrap);
    return to_return;
}

template<>
inline google::protobuf::Value value_cast<string, void>(const string& wrap) {
    google::protobuf::Value to_return;
    to_return.set_string_value(wrap);
    return to_return;
}

// TODO: more value casts for e.g. vectors and ints and embedded messages.

template<typename AnnotationType, typename Annotated>
AnnotationType get_annotation(const Annotated& annotated, const string& name) {
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
void set_annotation(Annotated* annotated, const string& name, const AnnotationType& annotation) {
    // Get ahold of the struct
    auto* annotation_struct = Annotation<Annotated>::get_mutable(annotated);
    
    // Set the key to the wrapped value
    (*annotation_struct->mutable_fields())[name] = value_cast(annotation);
}

template<typename Annotated>
void clear_annotation(Annotated* annotated, const string& name) {
    // Get ahold of the struct
    auto* annotation_struct = Annotation<Annotated>::get_mutable(annotated);
    // Clear out that field
    annotation_struct->mutable_fields()->clear(name);
}

}



#endif

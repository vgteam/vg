/// \file json2pb.cpp
/// Replicate the json2pb interface on top of Protobuf 3's built-in JSON code.

// Some code is derived from `json2pb`:

/*
 * Copyright (c) 2013 Pavel Shramov <shramov@mexmat.net>
 *
 * json2pb is free software; you can redistribute it and/or modify
 * it under the terms of the MIT license. See LICENSE for details.
 */

#include <json2pb.h>

#include <string>
#include <cassert>

#include <google/protobuf/util/json_util.h>
#include <jansson.h>

using google::protobuf::Message;

int json_dump_std_string(const char *buf, size_t size, void *data)
{
	std::string *s = (std::string *) data;
	s->append(buf, size);
	return 0;
}

void json2pb(Message &msg, FILE* fp) {
    // JSON is self-delimiting, but that's not useful if we don't know how to
    // parse it.
    
    // Protobuf has a JSON parser in it, but it's not exposed in a form that is
    // useful to us for pulling out single messages; it only wants to stream.

    // We continue using Jansson and parse to Jansson, then back to string, and
    // then to Protobuf.

    // Parse one JSON message from the file.
    json_t *root;
	json_error_t error;

	root = json_loadf(fp, JSON_DISABLE_EOF_CHECK, &error);

	if (!root)
		throw std::runtime_error(std::string("Load failed: ") + error.text);

	if (!json_is_object(root))
		throw std::runtime_error("Malformed JSON: not an object");
        
    // Dump back to string
    std::string r;
    json_dump_callback(root, json_dump_std_string, &r, 0);
   
    json_decref(root);
   
    // Parse to Protobuf.
    json2pb(msg, r);
}

void json2pb(Message &msg, const std::string& buf) {
    auto status = google::protobuf::util::JsonStringToMessage(buf, &msg);
    
    if (!status.ok()) {
        // This generally will happen if someone feeds in the wrong type of JSON.
        // TODO: It would be nice to be able to find the neme of the offending non-existent field.
        throw std::runtime_error("Could not deserialize " + msg.GetTypeName() + ": " + status.ToString());
    }
}

void json2pb(Message &msg, const char *buf, size_t size) {
	std::string buffer(buf, size);
    json2pb(msg, buffer);
}

std::string pb2json(const Message &msg) {
    // Set options to preserve field names and not camel case them
    google::protobuf::util::JsonPrintOptions opts;
    opts.preserve_proto_field_names = true;

	std::string buffer;
    auto status = google::protobuf::util::MessageToJsonString(msg, &buffer, opts);
    
    if (!status.ok()) {
        throw std::runtime_error("Could not serialize " + msg.GetTypeName() + ": " + status.ToString());
    }
    
    return buffer;
}

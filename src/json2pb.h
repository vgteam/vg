/// \file json2pb.h
/// Originally from `json2pb` by  Pavel Shramov.
/// The only part of that library that remains is the interface; the
/// implementation is now using the provided JSON code in Protobuf.
/// This file also contains a bunch of new vg-specific code.

/*
 * Copyright (c) 2013 Pavel Shramov <shramov@mexmat.net>
 *
 * json2pb is free software; you can redistribute it and/or modify
 * it under the terms of the MIT license. See LICENSE for details.
 */

#ifndef VG_JSON2PB_H_INCLUDED
#define VG_JSON2PB_H_INCLUDED

#include <string>
#include <cstdio>

namespace google {
namespace protobuf {
class Message;
}
}

void json2pb(google::protobuf::Message &msg, const std::string& buf);
void json2pb(google::protobuf::Message &msg, FILE* file);
void json2pb(google::protobuf::Message &msg, const char *buf, size_t size);
std::string pb2json(const google::protobuf::Message &msg);



#endif//VG_JSON2PB_H_INCLUDED

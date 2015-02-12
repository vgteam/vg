#ifndef JSON_HPP
#define JSON_HPP

#include <iostream>
#include <google/protobuf/message.h>
#include "pb2json.h"

void to_json(std::ostream& out, ::google::protobuf::Message& m);

#endif

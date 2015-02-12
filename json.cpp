#include "json.hpp"

void to_json(std::ostream& out, ::google::protobuf::Message& m) {
    char *json2 = pb2json(m);
    out << json2 << std::endl;
    free(json2);
}

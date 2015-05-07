#include <jansson.h>
#include <string>
#include <google/protobuf/descriptor.h>
#include <google/protobuf/message.h>
//using namespace google::protobuf;
char *pb2json(const google::protobuf::Message &msg);
char *pb2json(google::protobuf::Message *msg,const char *buf,int len);
static json_t *parse_msg(const google::protobuf::Message *msg);
static json_t *parse_repeated_field(const google::protobuf::Message *msg,const google::protobuf::Reflection * ref,const google::protobuf::FieldDescriptor *field);

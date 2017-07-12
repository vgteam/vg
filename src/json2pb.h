/*
 * Copyright (c) 2013 Pavel Shramov <shramov@mexmat.net>
 *
 * json2pb is free software; you can redistribute it and/or modify
 * it under the terms of the MIT license. See LICENSE for details.
 */

#ifndef VG_JSON2PB_H
#define VG_JSON2PB_H

#include <string>
#include <cstdio>
#include <functional>
#include <vector>
#include <stream.hpp>
#include <iostream>

namespace google {
namespace protobuf {
class Message;
}
}

void json2pb(google::protobuf::Message &msg, const char *buf, size_t size);
void json2pb(google::protobuf::Message &msg, FILE *fp);
std::string pb2json(const google::protobuf::Message &msg);

// It's handy to be able to stream in JSON via vg view for testing.
// This helper class takes this functionality from vg view -J and
// makes it more generic, so it can be used for other types than Graph.
// This class provides reading from JSON stream, as well as buffered
// conversion between input JSON stream and output protobuf or JSON stream.
template <class T>
class JSONStreamHelper
{
public:
    JSONStreamHelper(const std::string& file_name);
    ~JSONStreamHelper();
    // get a callback function that will read an object at a time from json stream.
    std::function<bool(T&)> get_read_fn();
    // read json stream (using above fn), and directly write to out in either
    // protobuf or json format. 
    int64_t write(std::ostream& out, bool json_out = false, int64_t buf_size = 1000);
private:
    FILE* _fp;
};


// Implementation of above:

template <class T>
inline JSONStreamHelper<T>::JSONStreamHelper(const std::string& file_name) {
    if (file_name == "-") {
        // Read standard input
        _fp = stdin;
    } else {
        // Open the file for reading
        _fp = fopen(file_name.c_str(), "r");
    }
    
    if (_fp == nullptr) {
        // We didn't manage to open the file. Complain and exit, before we try to dereference the pointer.
        std::cerr << "error:[JSONStreamHelper] could not open " << file_name << ": " << strerror(errno) << std::endl;
        exit(1);
    }
}

template <class T>
inline JSONStreamHelper<T>::~JSONStreamHelper() {
    if (_fp != stdin) {
        fclose(_fp);
    }
}
  
template <class T>
inline std::function<bool(T&)> JSONStreamHelper<T>::get_read_fn() {
    return [&](T& obj) -> bool {
        // zap protobuf object, since we want to overwrite and not append
        obj = T();
      
        // Check if the file ends now, and skip whitespace between records.
        char peeked;
        do {
            peeked = fgetc(this->_fp);
            if(peeked == EOF) {
                // File ended or otherwise errored. TODO: check for other
                // errors and complain.
                return false;
            }
        } while(isspace(peeked));
        // Put it back
        ungetc(peeked, this->_fp);
        
        // Now we know we have non-whitespace between here and EOF.
        // If it's not JSON, we want to die. So read it as JSON.
        json2pb(obj, this->_fp);
        
        // We read it successfully!
        return true;
    };
}

template<class T>
inline int64_t JSONStreamHelper<T>::write(std::ostream& out, bool json_out,
                                          int64_t buf_size) {    
    std::function<bool(T&)> reader = get_read_fn();        
    std::vector<T> buf;
    int64_t total = 0;
    bool good = true;
    std::function<T(uint64_t)> lambda = [&](uint64_t i) -> T {return buf[i];};
    while (good) {
        T obj;
        good = reader(obj);
        if (good) {
            buf.push_back(obj);
        }
        if (!good || buf.size() >= buf_size) {
            if (!json_out) {
                stream::write(out, buf.size(), lambda);
            } else {
                for (int i = 0; i < buf.size(); ++i) {
                    out << pb2json(buf[i]);
                }
            }
            total += buf.size();
            buf.clear();
        }
    }
    out.flush();
    return total;
}


#endif//VG_JSON2PB_H

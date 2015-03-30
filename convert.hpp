#ifndef __CONVERT_H
#define __CONVERT_H

#include <sstream>

// converts the string into the specified type, setting r to the converted
// value and returning true/false on success or failure
template<typename T>
bool convert(const std::string& s, T& r) {
    std::istringstream iss(s);
    iss >> r;
    return iss.eof() ? true : false;
}

template<typename T>
std::string convert(const T& r) {
    std::ostringstream iss;
    iss << r;
    return iss.str();
}

#endif

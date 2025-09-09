#include "log.hpp"
#include <iostream>
#include <vector>

namespace vg {

cerrWrapper basic_log(const std::string& context) {
    std::cerr << context << " ";
    return cerrWrapper(false);
}

cerrWrapper warning(const std::string& context) {
    std::cerr << "warning" << context << ": ";
    return cerrWrapper(false);
}

cerrWrapper fatal_error(const std::string& context) {
    std::cerr << "error" << context << ": ";
    return cerrWrapper(true);
}

}
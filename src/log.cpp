#include "log.hpp"
#include <iostream>
#include <vector>

namespace vg {

cerrWrapper basic_log(const std::string& context) {
    return cerrWrapper("[" + context + "] ", false);
}

cerrWrapper warning(const std::string& context) {
    return cerrWrapper("warning[" + context + "] ", false);
}

cerrWrapper fatal_error(const std::string& context) {
    return cerrWrapper("error[" + context + "] ", true);
}

}
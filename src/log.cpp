#include "log.hpp"
#include <iostream>
#include <vector>

namespace vg {

namespace logging {

cerrWrapper info(const std::string& context) {
    return cerrWrapper("[" + context + "] ", false);
}

cerrWrapper warn(const std::string& context) {
    return cerrWrapper("warning[" + context + "] ", false);
}

cerrWrapper error(const std::string& context) {
    return cerrWrapper("error[" + context + "] ", true);
}

}

}
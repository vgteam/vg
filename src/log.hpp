#include <string>
#include <iostream>

/** \file
 * log.hpp: defines a basic logging system for vg.
 * Shipped along with utility.hpp
 *
 * There are three options here:
 * - basic_log(): do nothing but append a context string to the front of
 *   the log message. This is useful for logging progress messages.
 *   This is a generic logging utility.
 * - warning(): emit a warning message with a given context
 * - fatal_error(): emit an error message with a given context
 * 
 * Errors/warnings are output to std::cerr as:
 * error<context>: message
 * warning<context>: message
 * where "context" is a string like "[vg inject]"
 * 
 * The error/warning functions exist to have a standardized error/warning format
 * across the vg subcommands. You don't *have* to use them elsewhere, but
 * it's a good idea to have consistent formatting.
 * 
 * These functions return a cerrWrapper object that behaves like a stream.
 * You can stream things to it with operator<<, and when the object goes out
 * of scope, it will exit the program with a failure code (in the case of
 * fatal_error) or do nothing (in the case of warning and basic_log).
 * 
 * Technically you can get the same effect by streaming to cerr, but this
 * way you don't have to remember to exit after logging an error,
 * and we can change format in a centralized place.
 */

namespace vg {

// https://stackoverflow.com/a/25615354/
class cerrWrapper {
private:
    std::ostream* str;
    bool exit_on_destruct;

public:
    cerrWrapper(bool exit_on_destruct) :
        exit_on_destruct(exit_on_destruct) {
        this->str = &std::cerr;
    }

    template <typename T>
    cerrWrapper& operator<<(const T& t) {
        if (str) {
            *str << t;
        }
        return *this;
    }

    cerrWrapper& operator<<(std::ostream& (*manip)(std::ostream&)) {
        *str << manip;
        return *this;
    }

    ~cerrWrapper() {
        if (exit_on_destruct) {
            exit(EXIT_FAILURE);
        }
    }
};

/// Log to cerr with a standard format
/// "context" is caller context, e.g. "[vg inject]"
cerrWrapper basic_log(const std::string& context);
/// Emit a warning with a standard format
/// "context" is caller context, e.g. "[vg inject]"
cerrWrapper warning(const std::string& context);
/// Error with a standard format
/// Once the cerrWrapper goes out of scope,
/// the program will exit with an error code.
/// "context" is caller context, e.g. "[vg inject]"
cerrWrapper fatal_error(const std::string& context);
}
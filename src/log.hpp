#include <string>
#include <iostream>

/** \file
 * log.hpp: defines a basic logging system for vg.
 * Shipped along with utility.hpp
 *
 * There are three options here:
 * - info(): do nothing but append a context string to the front of
 *   the log message. This is useful for logging progress messages.
 *   This is a generic logging utility.
 * - warn(): emit a warning message with a given context
 * - error(): emit an error message with a given context
 * 
 * Errors/warnings are output to std::cerr as:
 * error[<context>]: message
 * warning[<context>]: message
 * where "context" is a string like "vg inject"
 * 
 * The error/warning functions exist to have a standardized error/warning format
 * across the vg subcommands. You don't *have* to use them elsewhere, but
 * it's a good idea to have consistent formatting.
 * 
 * These functions return a cerrWrapper object that behaves like a stream.
 * You can stream things to it with operator<<, and when the object goes out
 * of scope, it will exit the program with a failure code (in the case of
 * error) or do nothing (in the case of warn() and info()).
 * 
 * Logging functions may either be accessed via logging::info/warn/error, or
 * via a Logger object that is initialized with a context string.
 */

namespace vg {

// https://stackoverflow.com/a/25615354/
class cerrWrapper {
private:
    bool exit_on_destruct;
    // How far to indent each time
    size_t indent_length;
    // Should we indent next time?
    bool at_start_of_line;

public:
    cerrWrapper(std::string prefix, bool exit_on_destruct) :
        exit_on_destruct(exit_on_destruct) {
        std::cerr << prefix;
        // Remember indentation level
        indent_length = prefix.size();
        at_start_of_line = false;
    }

    template <typename T>
    cerrWrapper& operator<<(const T& t) {
        std::cerr << t;
        return *this;
    }

    cerrWrapper& operator<<(std::ostream& (*manip)(std::ostream&)) {
        std::cerr << manip;
        return *this;
    }

    ~cerrWrapper() {
        if (exit_on_destruct) {
            exit(EXIT_FAILURE);
        }
    }
};

namespace logging {

/// Log to cerr with a standard format
/// "context" is caller context, e.g. "vg inject"
cerrWrapper info(const std::string& context);
/// Emit a warning with a standard format
/// "context" is caller context, e.g. "vg inject"
cerrWrapper warn(const std::string& context);
/// Error with a standard format
/// Once the cerrWrapper goes out of scope,
/// the program will exit with an error code.
/// "context" is caller context, e.g. "vg inject"
cerrWrapper error(const std::string& context);

}

/// Class to set up at the start of a file
/// when you want to generate a bunch of loggers
class Logger {
private:
    std::string context;

public:
    Logger(const std::string& context) : context(context) {}

    inline cerrWrapper info() const {
        return logging::info(context);
    }
    inline cerrWrapper warn() const {
        return logging::warn(context);
    }
    inline cerrWrapper error() const {
        return logging::error(context);
    }
};
}
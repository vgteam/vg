#include <string>
#include <iostream>

/** \file
 * log.hpp: defines a basic logging system for vg.
 * Shipped along with utility.hpp
 *
 * There are three options here:
 * - emit_with_indent(): assuming that the first line is already indented,
 *   emit a message with a given indent, wrapping lines at 80 columns.
 *   This is a generic logging utility.
 * - emit_warning(): emit a warning message with a given context, optionally
 *   pretty-printed with emit_with_indent().
 * - error_and_exit(): emit an error message with a given context, optionally
 *   pretty-printed with emit_with_indent(), and then exit with EXIT_FAILURE.
 * 
 * Errors/warnings are output to std::cerr as:
 * error<context>: message
 * warning<context>: message
 * where "context" is a string like "[vg inject]"
 * 
 * The error/warning functions exist to have a standardized error/warning format
 * across the vg subcommands. You don't *have* to use them elsewhere, but
 * it's a good idea to have consistent formatting.
 */

namespace vg {

void error_and_exit(const std::string& context, const std::string& message, bool pretty_print = true);
/// Warn the user with a standard format
/// "context" is the context in which the warning occurs, e.g. "[vg inject]"
/// If pretty_print is true, will reformat the message to fit within 80 columns
/// via emit_with_indent(): honors "\n" and "\t", adding extra newlines
void emit_warning(const std::string& context, const std::string& message, bool pretty_print = true);

/// Write a message to an output stream with indentation
/// Inserts "indent" spaces on the beginning of each line
/// Honors existing newlines, but otherwise adds newlines between words
/// such that everything gets printed with width no more than 80
/// Assumes that the first line is already indented
/// (because you put something there and you want others to line up)
void emit_with_indent(std::ostream& outstream, const std::string& message, std::size_t indent);

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
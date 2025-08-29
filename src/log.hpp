#include <string>

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

/// Error (and exit with an error code) with a standard format
/// If you're erroring because a file doesn't exist,
/// consider using require_exists() instead.
/// "context" is the context in which the error occurred, e.g. "[vg inject]"
/// If pretty_print is true, will reformat the message to fit within 80 columns
/// via emit_with_indent(): honors "\n" and "\t", adding extra newlines
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

}
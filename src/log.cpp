#include "log.hpp"
#include <iostream>
#include <vector>

namespace vg {

void error_and_exit(const std::string& context, const std::string& message, bool pretty_print) {
    #pragma omp critical (cerr)
    {
        std::cerr << "error" << context << ": ";
        if (pretty_print) {
            emit_with_indent(std::cerr, message, 7 + context.size());
        } else {
            std::cerr << message << std::endl;
        }
    }
    exit(EXIT_FAILURE);
}

void emit_warning(const std::string& context, const std::string& message, bool pretty_print) {
    #pragma omp critical (cerr)
    {
        std::cerr << "warning" << context << ": ";
        if (pretty_print) {
            emit_with_indent(std::cerr, message, 9 + context.size());
        } else {
            std::cerr << message << std::endl;
        }
    }
}

void emit_with_indent(std::ostream& outstream, const std::string& message, std::size_t indent) {
    std::vector<std::string> lines;
    // Break up by pre-existing newlines
    size_t start_pos = 0;
    size_t next_pos = 0;
    while ((next_pos = message.find('\n', start_pos)) != std::string::npos) {
        lines.push_back(message.substr(start_pos, next_pos - start_pos));
        start_pos = next_pos + 1;
    }
    // Add remaining bit
    lines.push_back(message.substr(start_pos, message.length()));
    // Output each line with the indent
    bool is_first_line = true;
    for (auto& line : lines) {
        // Find-and-replace tabs with 4 spaces
        // https://stackoverflow.com/a/5878802
        size_t pos = line.find("\t");
        while (pos != std::string::npos) {
            line.replace(pos, 1, "    ");
        }

        size_t cur_col = indent;
        start_pos = 0;
        next_pos = 0;
        if (!is_first_line) {
            // Assume only the first line is pre-indented
            outstream << std::string(indent, ' ');
        } else {
            is_first_line = false;
        }
        while ((next_pos = line.find(' ', start_pos)) != std::string::npos) {
            size_t word_length = next_pos - start_pos;
            if (word_length + cur_col > 80) {
                // If the next word would go past the end of the line, break
                outstream << "\n" << std::string(indent, ' ');
                cur_col = indent;
            }
            outstream << line.substr(start_pos, word_length) << " ";
            start_pos = next_pos + 1;
            cur_col += word_length + 1; // +1 for the space
        }
        if ((line.size() - start_pos) + cur_col > 80) {
            // If the last word would go past the end of the line, break
            outstream << "\n" << std::string(indent, ' ');
        }
        // Output the last bit of the line
        outstream << line.substr(start_pos) << "\n";
    }
    outstream.flush();
}

}
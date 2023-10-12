/**
 *\file
 * options.cpp: option parser system implementation
 */

#include "options.hpp"

namespace vg {
namespace subcommand {

void TickChainLink::reset_along_chain() {
    reset_along_chain_parent();
}

bool TickChainLink::tick_along_chain() {
    std::cerr << "Tick along chain at " << this << std::endl;
    return tick_along_chain_parent();
}

void TickChainLink::reset_chain() {
    reset_along_chain();
}

bool TickChainLink::tick_chain() {
    std::cerr << "Default tick chain at " << this << std::endl;
    return tick_along_chain();
}

bool TickChainLink::is_static() const {
    return true;
}

TickChainLink& TickChainLink::chain(TickChainLink& next) {
    std::cerr << "Chain " << this << " onto parent " << &next << std::endl;

    // Attach next to us
    next.reset_along_chain_parent = [&]() {
        this->reset_along_chain();
    };
    next.tick_along_chain_parent = [&]() {
        return this->tick_along_chain();
    };
    
    // And return it for a nice chain of chain calls.
    return next;
}

std::function<void(const std::function<void(void)>&)> TickChainLink::get_iterator() {
    return [&](const std::function<void(void)>& iteratee) {
        // Start
        reset_chain();
        
        do {
            // Run iteratee
            iteratee();
            // And tick the whole chain before running again
        } while(tick_chain());
    };
}

int get_option_id() {
    static int id = 10000;
    return id++;
}

template<>
const char* get_metavar<size_t>() {
    return "INT";
}

template<>
const char* get_metavar<int>() {
    return "INT";
}

template<>
const char* get_metavar<int8_t>() {
    return "INT";
}

template<>
const char* get_metavar<bool>() {
    return "BOOL";
}

template<>
const char* get_metavar<double>() {
    return "FLOAT";
}

template<>
const char* get_metavar<std::string>() {
    return "NAME";
}

BaseValuation::BaseValuation(const std::string& option) : option(option) {
    // Nothing to do!
}

const ValidatorFunction<double> double_is_positive = [](const double& d) {
    if (d <= 0) {
        throw std::domain_error("must be strictly positive");
    }
};

const ValidatorFunction<double> double_is_nonnegative = [](const double& d) {
    if (d < 0) {
        throw std::domain_error("cannot be negative");
    }
};

const ValidatorFunction<size_t> size_t_is_nonzero = [](const size_t& s) {
    if (s == 0) {
        throw std::domain_error("cannot be zero");
    }
};

const ValidatorFunction<int> int_is_nonnegative = [](const int& i) {
    if (i < 0) {
        throw std::domain_error("cannot be negative");
    }
};

TickChainLink& GroupedOptionGroup::chain(TickChainLink& next) {
    if (subgroups.empty()) {
        // Just chain through
        return TickChainLink::chain(next);
    } else {
        // We are already chained to first subgroup, so chain last subgroup to next.
        subgroups.back()->chain(next);
        return next;
    }
}

void GroupedOptionGroup::reset_chain() {
    if (subgroups.empty()) {
        TickChainLink::reset_chain();
    } else {
        // Delegate tick to the real end of the chain
        subgroups.back()->reset_chain();
    } 
}

bool GroupedOptionGroup::tick_chain() {
    std::cerr << "Grouped group tick chain at " << this << std::endl;
    if (!subgroups.empty()) {
        // Delegate tick to the real end of the chain
        return subgroups.back()->tick_chain();
    }
    return false;
}

bool GroupedOptionGroup::parse(int option_id, const char* optarg) {
    for (auto& group : subgroups) {
        if (group->parse(option_id, optarg)) {
            // If any of our groups wants this option, we do too.
            return true;
        }
    }
    return false;
}

bool GroupedOptionGroup::preset(const BaseValuation& entry) {
    for (auto& group : subgroups) {
        if (group->preset(entry)) {
            // If any of our groups wants this option, we do too.
            return true;
        }
    }
    return false;
}

bool GroupedOptionGroup::set(const BaseValuation& entry) {
    for (auto& group : subgroups) {
        if (group->set(entry)) {
            // If any of our groups wants this option, we do too.
            return true;
        }
    }
    return false;
}

bool GroupedOptionGroup::query(BaseValuation& entry) const {
    for (auto& group : subgroups) {
        if (group->query(entry)) {
            // If any of our groups wants this option, we do too.
            return true;
        }
    }
    return false;
}

void GroupedOptionGroup::print_options(ostream& out, bool slug) const {
    for (auto& group : subgroups) {
        // Print options from all groups in order
        group->print_options(out, slug);
    }
}

std::vector<std::pair<std::string, std::string>> GroupedOptionGroup::get_help() const {
    std::vector<std::pair<std::string, std::string>> helps;
    for (auto& group : subgroups) {
        // Get helps from all subgroups
        auto subgroup_helps = group->get_help();
        // And put each collection on the end
        std::copy(subgroup_helps.begin(), subgroup_helps.end(), std::back_inserter(helps));
    }
    return helps;
}

void GroupedOptionGroup::make_long_options(std::vector<struct option>& dest) const {
    for (auto& group : subgroups) {
        group->make_long_options(dest);
    }
}

void GroupedOptionGroup::make_short_options(std::string& dest) const {
    for (auto& group : subgroups) {
        group->make_short_options(dest);
    }
}

void print_table(const std::vector<std::pair<std::string, std::string>>& rows, ostream& out) {
    // Work out the max length of anything in the first column
    size_t max_length = 0;
    for (auto& r : rows) {
        max_length = std::max(max_length, r.first.size());
    }
    for (auto& r : rows) {
        if (r.first.empty()) {
            // It's a heading
            out << r.second << std::endl;
        } else {
            // Print leading indent
            out << "  ";
            // Print column 1
            out << r.first;
            for (size_t i = 0; i < max_length - r.first.size(); i++) {
                // Print padding to make all items the max length
                out << " ";
            }
            // Print separator
            out << " ";
            // Print column 2
            out << r.second << std::endl;
        }
    }
}

}
}

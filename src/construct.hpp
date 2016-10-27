#ifndef VG_CONSRUCT_HPP
#define VG_CONSTRUCT_HPP

/*
 * construct.hpp: contains the subcommand implementation for `vg construct`, and
 * the actual vg-from-vcf construction code.
 */

void help_construct(char** argv);
int main_construct(int argc, char** argv);

// Actual constructing-things code is part of vg::VG, with declarations in
// vg.hpp.

#endif

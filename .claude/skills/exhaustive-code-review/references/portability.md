# Portability lens

Ask only: does this build and behave the same everywhere vg is built?

vg builds on Linux x86-64 and ARM64 with GCC, and on macOS with Apple Clang
(ARM64). The Makefile and `source_me.sh` handle the differences; the CI files
under `.github/` and `.gitlab-ci.yml` show which configurations are tested.

- **Compiler differences.** GCC-only extensions, a missing `#include` that one
  standard library pulls in by accident (`<cstdint>`, `<algorithm>`,
  `<functional>`), `std::` functions not yet in Apple's library, and C++
  features newer than the standard vg sets in `CXX_STANDARD` in the Makefile.
- **Sizes.** `long` is 64 bits on both platforms here, but code should use
  `int64_t`/`size_t` for anything that can exceed 32 bits. Node IDs are
  `nid_t` (64-bit).
- **Endianness and layout.** Code that writes structs to disk or reads them
  back with `memcpy` assumes layout and byte order. vg's formats go through
  Protobuf or libraries' own serializers; check that new on-disk data does too.
- **Platform calls.** `/proc`, `sysconf`, `sched_getaffinity`, `malloc_usable_size`
  and Linux-only headers need a macOS path or a guard. File-system case
  sensitivity differs on macOS.
- **Build system.** New source files are picked up automatically. A new
  dependency needs to work in the Makefile's dependency rules, on macOS with
  Homebrew, and in the Docker build.
- **Randomness and hashing.** Output that depends on `std::hash` or
  unordered-container iteration order can differ between standard libraries.
  Tests that compare such output break on one platform only.

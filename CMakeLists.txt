cmake_minimum_required(VERSION 3.22)

project(vg
    LANGUAGES CXX C
    DESCRIPTION "Variation Graph Tools"
)

# ────────────────────────────────────────────────────────────────
# Output directories (mirror the Makefile's bin/ lib/ layout)
# ────────────────────────────────────────────────────────────────
# ────────────────────────────────────────────────────────────────
# Submodule check — fail early with a clear message if deps/ is empty
# ────────────────────────────────────────────────────────────────
if(NOT EXISTS ${CMAKE_SOURCE_DIR}/deps/htslib/configure.ac)
    message(FATAL_ERROR
        "Git submodules appear to be uninitialized.\n"
        "Please run:\n"
        "  git submodule update --init --recursive\n"
        "before running CMake."
    )
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# Internal staging directories used by dep builds
set(VG_LIB_DIR     ${CMAKE_BINARY_DIR}/lib)
set(VG_INC_DIR     ${CMAKE_BINARY_DIR}/include)
set(VG_BIN_DIR     ${CMAKE_BINARY_DIR}/bin)
set(VG_PKG_DIR     ${CMAKE_BINARY_DIR}/lib/pkgconfig)

file(MAKE_DIRECTORY ${VG_LIB_DIR} ${VG_INC_DIR} ${VG_BIN_DIR} ${VG_PKG_DIR})

# ────────────────────────────────────────────────────────────────
# C++ standard — C++20 everywhere (gbwtgraph uses inline variables)
# ────────────────────────────────────────────────────────────────
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

#-───────────────────────────────────────────────────────────────
set(CMAKE_POLICY_VERSION_MINIMUM 3.5)

# ────────────────────────────────────────────────────────────────
# Build variant options  (replicate Makefile's asan/jemalloc/mimalloc vars)
# ────────────────────────────────────────────────────────────────
option(VG_ASAN        "Build with AddressSanitizer (disables jemalloc)" OFF)
option(VG_USE_JEMALLOC "Use jemalloc allocator"  ON)
option(VG_USE_MIMALLOC "Use mimalloc allocator"  OFF)
set(RUN_DOXYGEN OFF)

# Mutual-exclusion validation
if(VG_USE_JEMALLOC AND VG_USE_MIMALLOC)
    message(FATAL_ERROR "VG_USE_JEMALLOC and VG_USE_MIMALLOC cannot both be ON")
endif()
if(VG_ASAN AND VG_USE_JEMALLOC)
    message(WARNING "ASAN requested — disabling jemalloc (set VG_USE_JEMALLOC=OFF to suppress this warning)")
    set(VG_USE_JEMALLOC OFF CACHE BOOL "" FORCE)
endif()

# macOS never uses jemalloc (mirrors Makefile: ifeq $(uname -s),Darwin → jemalloc=off)
if(APPLE AND VG_USE_JEMALLOC)
    message(STATUS "macOS detected — disabling jemalloc (not supported on this platform)")
    set(VG_USE_JEMALLOC OFF CACHE BOOL "" FORCE)
endif()

# ────────────────────────────────────────────────────────────────
# Platform and compiler detection
# ────────────────────────────────────────────────────────────────
if(APPLE)
    message(STATUS "OS is Mac")
else()
    message(STATUS "OS is Linux")
endif()

# Detect Clang vs GCC (needed for OpenMP and libstdc++ handling on Mac)
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(VG_COMPILER_IS_CLANG TRUE)
    message(STATUS "Compiler ${CMAKE_CXX_COMPILER} is Clang")
else()
    set(VG_COMPILER_IS_CLANG FALSE)
    message(STATUS "Compiler ${CMAKE_CXX_COMPILER} is assumed to be GCC")
endif()

# ────────────────────────────────────────────────────────────────
# Base compile options  (replicate Makefile CXXFLAGS baseline)
# ────────────────────────────────────────────────────────────────
add_compile_options(
    -O3
    -Werror=return-type
    -ggdb
    -g
)

if(VG_ASAN)
    add_compile_options(-fsanitize=address -fno-omit-frame-pointer)
    add_link_options(-fsanitize=address)
endif()

# Architecture-specific flags
if(APPLE)
    if(CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64")
        add_compile_options(-march=native)
    endif()
    # -march=native is not recognized by Apple Clang on ARM; skip it there
else()
    # Linux: build for SSE4.2 for vaguely portable binaries
    if(CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64")
        add_compile_options(-msse4.2)
    endif()
    # Linux: allow backtrace access to all symbols
    add_link_options(-rdynamic)
endif()

# ────────────────────────────────────────────────────────────────
# OpenMP
# ────────────────────────────────────────────────────────────────
if(APPLE AND VG_COMPILER_IS_CLANG)
    # Apple Clang doesn't bundle libomp; use -Xpreprocessor -fopenmp
    # and find libomp from Homebrew or Macports.
    find_package(OpenMP REQUIRED)
    # CMake 3.16+ sets OpenMP_CXX_FLAGS correctly for AppleClang,
    # but we verify the omp header is reachable.
    add_compile_options(-Xpreprocessor -fopenmp)
    # OpenMP_CXX_LIBRARIES is set by find_package
    set(VG_OMP_LINK_LIBS ${OpenMP_CXX_LIBRARIES})
else()
    find_package(OpenMP REQUIRED)
    add_compile_options(-fopenmp)
    set(VG_OMP_LINK_LIBS ${OpenMP_CXX_LIBRARIES})
endif()

# ────────────────────────────────────────────────────────────────
# pkg-config system dependencies
# ────────────────────────────────────────────────────────────────
find_package(PkgConfig REQUIRED)

# Dynamic system deps (replicate PKG_CONFIG_DEPS := cairo libzstd)
pkg_check_modules(CAIRO   REQUIRED IMPORTED_TARGET cairo)
pkg_check_modules(ZSTD    REQUIRED IMPORTED_TARGET libzstd)

# Static system deps (replicate PKG_CONFIG_STATIC_DEPS := protobuf jansson)
find_package(Protobuf CONFIG REQUIRED)
pkg_check_modules(JANSSON REQUIRED IMPORTED_TARGET jansson)

# ────────────────────────────────────────────────────────────────
# Boost Program Options  (Makefile: -lboost_program_options[-mt])
# ────────────────────────────────────────────────────────────────
find_package(Boost REQUIRED COMPONENTS program_options)

# Python3 interpreter — needed by vcflib's pybind11 and the vgmanmd.py manpage generator
find_package(Python3 COMPONENTS Interpreter)

# ────────────────────────────────────────────────────────────────
# System library flags — collected into a CMake INTERFACE target
# so any target that needs them can link against vg_system_libs.
# ────────────────────────────────────────────────────────────────
add_library(vg_system_libs INTERFACE)

target_link_libraries(vg_system_libs INTERFACE
    PkgConfig::CAIRO
    PkgConfig::ZSTD
    protobuf::libprotobuf
    PkgConfig::JANSSON
    Boost::program_options
    ${VG_OMP_LINK_LIBS}
    pthread
    ncurses
    crypto
    dl
    lzma
    bz2
)

if(APPLE)
    # Mac: nothing extra (elfutils is Linux-only; no -latomic needed)
else()
    # Linux: elfutils libraries (built as dep in Phase 5) + atomic
    target_link_libraries(vg_system_libs INTERFACE
        atomic
        # elfutils targets added in Phase 5:
        #   dwfl dw dwelf ebl elf
    )
endif()

# Static section libraries (replicate LD_STATIC_LIB_FLAGS)
# These will be linked statically in the final executable link command.
# Listed here as a variable; applied in Phase 8 with appropriate
# -Wl,-Bstatic / -Wl,-Bdynamic wrappers on Linux.
set(VG_STATIC_LIBS
    # libvgio hts deflate z bz2 lzma — added when dep targets are built
)

# ────────────────────────────────────────────────────────────────
# Homebrew / Macports prefix hints (passed to CMake deps as CMAKE_PREFIX_PATH)
# ────────────────────────────────────────────────────────────────
if(APPLE)
    execute_process(
        COMMAND brew --prefix
        OUTPUT_VARIABLE HOMEBREW_PREFIX
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
    )
    if(HOMEBREW_PREFIX)
        message(STATUS "Homebrew prefix: ${HOMEBREW_PREFIX}")
        list(APPEND CMAKE_PREFIX_PATH ${HOMEBREW_PREFIX})
        include_directories(SYSTEM ${HOMEBREW_PREFIX}/include)
        link_directories(${HOMEBREW_PREFIX}/lib)

        # Bison from Homebrew (replaces Apple's ancient Bison)
        set(BISON_EXECUTABLE ${HOMEBREW_PREFIX}/opt/bison/bin/bison)
    endif()
    # Macports
    if(EXISTS /opt/local/lib)
        link_directories(/opt/local/lib)
    endif()
    if(EXISTS /usr/local/lib)
        link_directories(/usr/local/lib)
    endif()
endif()

# OMP_PREFIXES used by libvgio and vcflib's CMake to find OpenMP
set(VG_OMP_PREFIXES "/")
if(APPLE AND HOMEBREW_PREFIX)
    if(EXISTS ${HOMEBREW_PREFIX}/include/omp.h)
        list(APPEND VG_OMP_PREFIXES ${HOMEBREW_PREFIX})
    elseif(EXISTS ${HOMEBREW_PREFIX}/opt/libomp/include)
        list(APPEND VG_OMP_PREFIXES ${HOMEBREW_PREFIX}/opt/libomp)
    elseif(EXISTS /opt/local/lib/libomp)
        list(APPEND VG_OMP_PREFIXES /opt/local)
    endif()
endif()

# Escape semicolons so the list can be safely embedded inside a quoted CMAKE_ARGS string
# (e.g. "-DCMAKE_PREFIX_PATH=${CMAKE_BINARY_DIR}\;/usr\;${VG_OMP_PREFIXES_STR}")
string(REPLACE ";" "\\;" VG_OMP_PREFIXES_STR "${VG_OMP_PREFIXES}")

# ────────────────────────────────────────────────────────────────
# Shared-library suffix (used in dep build rules; CMake usually
# handles this automatically via CMAKE_SHARED_LIBRARY_SUFFIX)
# ────────────────────────────────────────────────────────────────
if(APPLE)
    set(VG_SHARED_SUFFIX dylib)
else()
    set(VG_SHARED_SUFFIX so)
endif()

# Static linking wrappers (Linux only; empty on Mac)
if(APPLE)
    set(VG_START_STATIC "")
    set(VG_END_STATIC "")
else()
    set(VG_START_STATIC -Wl,-Bstatic)
    set(VG_END_STATIC   -Wl,-Bdynamic)
endif()

# ────────────────────────────────────────────────────────────────
# Position-independent code — required for the shared libvg
# ────────────────────────────────────────────────────────────────
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# ────────────────────────────────────────────────────────────────
# Phases 3–11 are added via subdirectory includes below.
# Each include() corresponds to one conversion phase.
# Uncomment each line as that phase is implemented.
# ────────────────────────────────────────────────────────────────
include(cmake/deps_header_only.cmake)        # Phase 3
include(cmake/version_headers.cmake)         # Phase 10 (run early: generates src/ headers before compilation)
include(cmake/deps_cmake_native.cmake)       # Phase 4
include(cmake/deps_external_project.cmake)   # Phase 5
include(cmake/deps_xg.cmake)                 # Phase 6
include(cmake/libvg.cmake)                   # Phase 7
include(cmake/vg_executable.cmake)           # Phase 8
include(cmake/vg_tests.cmake)                # Phase 9
include(cmake/utility_targets.cmake)         # Phase 11

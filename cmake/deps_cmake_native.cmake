# Phase 4: CMake-based bundled dependencies
#
# Uses ExternalProject_Add for all of them (not add_subdirectory).
# Reasons:
#   - Prevents their install targets from polluting the parent CMake install
#   - Keeps their cache entries out of the parent build tree
#   - Lets us control the install prefix (VG_INC_DIR / VG_LIB_DIR)
#
# Deps covered here (standalone — no dependency on Phase 5 autoconf outputs):
#   kff-cpp-api, mimalloc (conditional), libhandlegraph, raptor
#
# Deps covered at the BOTTOM of cmake/deps_external_project.cmake (Phase 5)
# because they depend on htslib or tabixpp:
#   libvgio  (needs htslib)
#   vcflib   (needs htslib + tabixpp)
#   libbdsg  (needs libhandlegraph + sdsl)

include(ExternalProject)

set(DEPS_DIR ${CMAKE_SOURCE_DIR}/deps)

# Common cmake args propagated to all sub-builds
set(VG_DEP_CMAKE_ARGS
    -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
    -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
    -DCMAKE_CXX_STANDARD=17
    # Install directly into the build tree's lib/include dirs (matches VG_LIB_DIR/VG_INC_DIR)
    -DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}
    -DCMAKE_INSTALL_LIBDIR=lib
    -DCMAKE_POSITION_INDEPENDENT_CODE=ON
    -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
)

# ── kff-cpp-api ───────────────────────────────────────────────────────────
# Makefile: cmake -DCMAKE_CXX_FLAGS="-fPIC -Wall -Ofast -g $(CXXFLAGS)" ..
#            && make && cp kff_io.hpp <inc> && mv libkff.a <lib>
add_subdirectory(${DEPS_DIR}/kff-cpp-api ${CMAKE_BINARY_DIR}/build/kff EXCLUDE_FROM_ALL)
add_library(dep_kff STATIC IMPORTED GLOBAL)
set_target_properties(dep_kff PROPERTIES
    IMPORTED_LOCATION ${VG_LIB_DIR}/libkff.a
)
target_link_libraries(dep_kff INTERFACE kff)

# ── mimalloc (optional) ───────────────────────────────────────────────────
# Makefile: cmake ... && make && cp include/* <inc>/ && cp mimalloc.o <lib>/
if(VG_USE_MIMALLOC)
    ExternalProject_Add(mimalloc_ep
        SOURCE_DIR ${DEPS_DIR}/mimalloc
        BINARY_DIR ${CMAKE_BINARY_DIR}/build/mimalloc
        CMAKE_ARGS
            ${VG_DEP_CMAKE_ARGS}
            -DCMAKE_C_FLAGS=${CMAKE_C_FLAGS}
            -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}
        BUILD_COMMAND     ${CMAKE_MAKE_PROGRAM}
        INSTALL_COMMAND
            ${CMAKE_COMMAND} -E copy_directory
                ${DEPS_DIR}/mimalloc/include ${VG_INC_DIR}
            COMMAND ${CMAKE_COMMAND} -E copy
                <BINARY_DIR>/mimalloc.o ${VG_LIB_DIR}/mimalloc.o
        BUILD_BYPRODUCTS ${VG_LIB_DIR}/mimalloc.o
    )
    # mimalloc is linked as an object file, not a library; exposed as a var
    set(VG_MIMALLOC_OBJ ${VG_LIB_DIR}/mimalloc.o)
    add_custom_target(dep_mimalloc DEPENDS mimalloc_ep)
endif()

# ── libhandlegraph ────────────────────────────────────────────────────────
# Makefile: cmake -DCMAKE_INSTALL_PREFIX=$(CWD) -DCMAKE_INSTALL_LIBDIR=lib ..
#            && make && make install
add_subdirectory(${DEPS_DIR}/libhandlegraph ${CMAKE_BINARY_DIR}/build/libhandlegraph EXCLUDE_FROM_ALL)
# Stage handlegraph headers into the shared include dir so Makefile-based deps can find them
add_custom_target(handlegraph_stage_headers
    COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${DEPS_DIR}/libhandlegraph/src/include/handlegraph
        ${VG_INC_DIR}/handlegraph
)
#add_dependencies(handlegraph_stage_headers libhandlegraph)

add_library(dep_libhandlegraph STATIC IMPORTED GLOBAL)
set_target_properties(dep_libhandlegraph PROPERTIES
    IMPORTED_LOCATION             ${VG_LIB_DIR}/libhandlegraph.a
    INTERFACE_INCLUDE_DIRECTORIES ${VG_INC_DIR}
)
target_link_libraries(dep_libhandlegraph INTERFACE handlegraph_stage_headers)

# ── raptor ────────────────────────────────────────────────────────────────
# Makefile recipe (complex — must regenerate turtle lexer from Bison):
#   cd raptor/build
#   rm CMakeCache.txt CMakeFiles ...
#   cmake -DCMAKE_POLICY_VERSION_MINIMUM=3.5 ...
#   rm -f src/turtle_parser.c src/turtle_lexer.c
#   make turtle_lexer_tgt
#   make -f src/CMakeFiles/raptor2.dir/build.make src/turtle_lexer.c
#   sed -i.bak '/yycleanup/d' src/turtle_lexer.c
#   make
#   cp src/libraptor2.a <lib>/
#   cp build/utils/rapper <bin>/
#   cp build/src/*.h <inc>/raptor2/

# We need Bison; on Mac use Homebrew's version (set in Phase 2 scaffold)
find_program(BISON_EXECUTABLE bison
    HINTS ${HOMEBREW_PREFIX}/opt/bison/bin /usr/local/bin /usr/bin
    REQUIRED
)
message(STATUS "Bison: ${BISON_EXECUTABLE}")

# Build command script: configure + lexer hack + build
ExternalProject_Add(raptor_ep
    SOURCE_DIR ${DEPS_DIR}/raptor
    BINARY_DIR ${DEPS_DIR}/raptor/build      # Makefile uses in-source build dir
    # No configure step; we do it in BUILD_COMMAND so we can also do the lexer hack
    CONFIGURE_COMMAND ""
    BUILD_COMMAND
        ${CMAKE_COMMAND} -E remove -f
            <BINARY_DIR>/CMakeCache.txt
            <BINARY_DIR>/src/turtle_parser.c
            <BINARY_DIR>/src/turtle_lexer.c
        COMMAND
            ${CMAKE_COMMAND}
                -DCMAKE_POLICY_VERSION_MINIMUM=3.5
                -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                -DCMAKE_C_FLAGS=-fPIC\ ${CMAKE_C_FLAGS}
                -DCMAKE_CXX_FLAGS=-fPIC\ ${CMAKE_CXX_FLAGS}
                -DBISON_EXECUTABLE=${BISON_EXECUTABLE}
                -B <BINARY_DIR>
                -S <SOURCE_DIR>
        COMMAND ${CMAKE_MAKE_PROGRAM} -C <BINARY_DIR> turtle_lexer_tgt
        COMMAND ${CMAKE_MAKE_PROGRAM} -C <BINARY_DIR>
            -f src/CMakeFiles/raptor2.dir/build.make
            src/turtle_lexer.c
        COMMAND
            # Remove yycleanup — causes build errors (Makefile uses sed -i.bak)
            ${CMAKE_COMMAND} -P
            ${CMAKE_SOURCE_DIR}/cmake/raptor_patch_lexer.cmake
        COMMAND ${CMAKE_MAKE_PROGRAM} -C <BINARY_DIR>
    INSTALL_COMMAND
        ${CMAKE_COMMAND} -E make_directory ${VG_INC_DIR}/raptor2
        COMMAND ${CMAKE_COMMAND} -E copy
            <BINARY_DIR>/src/libraptor2.a ${VG_LIB_DIR}/libraptor2.a
        COMMAND ${CMAKE_COMMAND} -E copy
            <BINARY_DIR>/utils/rapper     ${VG_BIN_DIR}/rapper
        COMMAND ${CMAKE_COMMAND} -E copy_directory
            <BINARY_DIR>/src             ${VG_INC_DIR}/raptor2
            # (copies *.h; CMake copy_directory copies all files)
    BUILD_BYPRODUCTS
        ${VG_LIB_DIR}/libraptor2.a
        ${VG_INC_DIR}/raptor2/raptor2.h
        ${VG_BIN_DIR}/rapper
    DEPENDS ""  # No cross-deps; standalone
)

add_library(dep_raptor STATIC IMPORTED GLOBAL)
set_target_properties(dep_raptor PROPERTIES
    IMPORTED_LOCATION             ${VG_LIB_DIR}/libraptor2.a
    INTERFACE_INCLUDE_DIRECTORIES ${VG_INC_DIR}
)
add_dependencies(dep_raptor raptor_ep)

# ── raptor_patch_lexer.cmake helper ───────────────────────────────────────
# Writes a small cmake -P script that replicates the sed command.
# We generate it at configure time so it can embed the correct binary path.
file(WRITE ${CMAKE_SOURCE_DIR}/cmake/raptor_patch_lexer.cmake
"# Generated by cmake/deps_cmake_native.cmake
# Replicates: sed -i.bak '/yycleanup/d' raptor/build/src/turtle_lexer.c
set(LEXER_FILE \"${DEPS_DIR}/raptor/build/src/turtle_lexer.c\")
if(EXISTS \${LEXER_FILE})
    file(READ \${LEXER_FILE} _content)
    string(REGEX REPLACE \"[^\\n]*yycleanup[^\\n]*\\n\" \"\" _patched \"\${_content}\")
    file(WRITE \${LEXER_FILE} \"\${_patched}\")
    message(STATUS \"Patched raptor turtle_lexer.c (removed yycleanup)\")
endif()
")

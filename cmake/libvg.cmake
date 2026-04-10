# Phase 7: libvg — static and shared library targets
#
# Mirrors the Makefile's:
#   $(LIB_DIR)/libvg.a     — static archive (no linking, just ar)
#   $(LIB_DIR)/libvg.so/dylib — linked shared library
#
# Sources:
#   src/*.cpp              (minus main.cpp)  → OBJ
#   src/algorithms/*.cpp                     → ALGORITHMS_OBJ
#   src/io/*.cpp                             → IO_OBJ
#   deps/sha1/sha1.cpp                       → DEP_OBJ (via dep_sha1 OBJECT lib)
#   deps/progress_bar/progress_bar.cpp       → DEP_OBJ (via dep_progress_bar OBJECT lib)

set(VG_SRC_DIR  ${CMAKE_SOURCE_DIR}/src)

# ── Source file globs ──────────────────────────────────────────────────────
file(GLOB VG_SRCS         ${VG_SRC_DIR}/*.cpp)
file(GLOB ALG_SRCS        ${VG_SRC_DIR}/algorithms/*.cpp)
file(GLOB IO_SRCS         ${VG_SRC_DIR}/io/*.cpp)

# Exclude main.cpp — it goes only in the vg executable (Phase 8)
list(REMOVE_ITEM VG_SRCS ${VG_SRC_DIR}/main.cpp)

# ── Include directories (mirrors Makefile's INCLUDE_FLAGS) ────────────────
# -I$(CWD)/$(INC_DIR) -I. -I$(CWD)/$(SRC_DIR)
# -I$(CWD)/$(UNITTEST_SRC_DIR) -I$(CWD)/$(UNITTEST_SUPPORT_SRC_DIR)
# -I$(CWD)/$(SUBCOMMAND_SRC_DIR) -I$(CWD)/$(INC_DIR)/dynamic
# + pkg-config INCLUDE_FLAGS from system deps

set(VG_INCLUDE_DIRS
    ${VG_INC_DIR}                                # staging include (dep headers)
    ${VG_INC_DIR}/dynamic                         # DYNAMIC's own include subdir
    ${CMAKE_SOURCE_DIR}                           # -I. (project root)
    ${VG_SRC_DIR}                                 # -I src/
    ${VG_SRC_DIR}/unittest                        # -I src/unittest
    ${VG_SRC_DIR}/unittest/support                # -I src/unittest/support
    ${VG_SRC_DIR}/subcommand                      # -I src/subcommand
)

# ── Static library: libvg.a ───────────────────────────────────────────────
# The Makefile's libvg.a is a plain archive (no linking step).
# In CMake, add_library(STATIC) produces exactly that.

add_library(libvg STATIC
    ${VG_SRCS}
    ${ALG_SRCS}
    ${IO_SRCS}
    # DEP_OBJ: sha1 and progress_bar OBJECT libs are added via generator expressions
    $<TARGET_OBJECTS:dep_sha1>
    $<TARGET_OBJECTS:dep_progress_bar>
)

set_target_properties(libvg PROPERTIES
    OUTPUT_NAME vg
    POSITION_INDEPENDENT_CODE ON
)

target_include_directories(libvg PUBLIC ${VG_INCLUDE_DIRS})

# vg source uses POSIX-era system headers; ensure C++17 is set
target_compile_features(libvg PUBLIC cxx_std_17)

# Pass the git/env version headers through (generated in Phase 10;
# add a dependency so CMake re-runs if vg_git_version.hpp changes)
target_sources(libvg PRIVATE ${VG_SRC_DIR}/version.cpp)
# (version.cpp is already in VG_SRCS from the glob, but listing explicitly
#  makes the dependency visible for documentation purposes)

# Link against all dep IMPORTED targets — these express include paths
# (via INTERFACE_INCLUDE_DIRECTORIES) but no actual linking happens yet;
# the linker flags come in when libvg itself is linked into an executable.
target_link_libraries(libvg PUBLIC
    # Phase 3: header-only
    dep_headers_all

    # Phase 4: CMake-native deps
    dep_libhandlegraph
    dep_kff
    dep_raptor

    # Phase 5: non-CMake deps
    dep_libdeflate
    dep_ssw
    dep_sdsl
    dep_htslib
    dep_snappy
    dep_sparsehash
    dep_gcsa2
    dep_gbwt
    dep_gbwtgraph
    dep_tabixpp
    dep_gssw
    dep_sonlib
    dep_pinchesandcacti
    dep_3edgeconnected
    dep_libvcfh
    dep_fml
    dep_structures
    dep_sublinearLS
    dep_vcflib
    dep_wfa2
    dep_libvgio
    dep_libbdsg

    # Phase 6: xg
    dep_xg

    # System deps
    vg_system_libs
)

# On Linux, add elfutils dep targets (only defined on non-Apple)
if(NOT APPLE)
    target_link_libraries(libvg PUBLIC dep_elf dep_ebl dep_dwfl dep_dwelf dep_dw)
endif()

# Ensure all ExternalProject targets are built before libvg sources compile
add_dependencies(libvg
    dep_libhandlegraph kff raptor_ep
    libdeflate_ep ssw_ep dep_sdsl htslib_ep snappy_ep sparsehash_ep
    gcsa2_ep gbwt_ep gbwtgraph_ep tabixpp_ep gssw_ep sonlib_ep
    pinchescacti_ep libvcfh_ep fml_ep structures_ep sublinearls_ep
    vcflib_ep dep_libvgio dep_libbdsg xg_ep
)
if(NOT APPLE)
    add_dependencies(libvg elfutils_ep)
endif()
if(VG_USE_JEMALLOC)
    add_dependencies(libvg jemalloc_ep)
endif()

# ── Shared library: libvg.so / libvg.dylib ───────────────────────────────
# The Makefile links the shared lib with -shared and all dep libraries.
# In CMake, SHARED libraries are linked (not just archived).

add_library(libvg_shared SHARED
    ${VG_SRCS}
    ${ALG_SRCS}
    ${IO_SRCS}
    $<TARGET_OBJECTS:dep_sha1>
    $<TARGET_OBJECTS:dep_progress_bar>
)

set_target_properties(libvg_shared PROPERTIES
    OUTPUT_NAME     vg
    POSITION_INDEPENDENT_CODE ON
    # Avoid name clash with the STATIC target's output
    LIBRARY_OUTPUT_DIRECTORY ${VG_LIB_DIR}
)

target_include_directories(libvg_shared PUBLIC ${VG_INCLUDE_DIRS})
target_compile_features(libvg_shared PUBLIC cxx_std_17)

# Shared lib needs the same deps, plus it actually links them at build time
target_link_libraries(libvg_shared PRIVATE
    dep_headers_all
    dep_libhandlegraph dep_kff dep_raptor
    dep_libdeflate dep_ssw dep_sdsl
    dep_htslib dep_snappy dep_sparsehash
    dep_gcsa2 dep_gbwt dep_gbwtgraph dep_tabixpp dep_gssw
    dep_sonlib dep_pinchesandcacti dep_3edgeconnected dep_libvcfh
    dep_fml dep_structures dep_sublinearLS
    dep_vcflib dep_wfa2 dep_libbdsg
    dep_xg
    vg_system_libs
)
if(NOT APPLE)
    target_link_libraries(libvg_shared PRIVATE dep_elf dep_ebl dep_dwfl dep_dwelf dep_dw)
endif()

add_dependencies(libvg_shared
    dep_libhandlegraph kff raptor_ep
    libdeflate_ep ssw_ep dep_sdsl htslib_ep snappy_ep sparsehash_ep
    gcsa2_ep gbwt_ep gbwtgraph_ep tabixpp_ep gssw_ep sonlib_ep
    pinchescacti_ep libvcfh_ep fml_ep structures_ep sublinearls_ep
    vcflib_ep dep_libvgio dep_libbdsg xg_ep
)
if(NOT APPLE)
    add_dependencies(libvg_shared elfutils_ep)
endif()

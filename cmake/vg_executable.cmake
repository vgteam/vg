# Phase 8: vg executable
#
# Mirrors the Makefile's final link step:
#
#   $(CXX) ... -o bin/vg
#     [mimalloc.o]                      # PRE_LINK_DEPS (if mimalloc=on)
#     obj/main.o                        # src/main.cpp
#     obj/unittest/*.o                  # src/unittest/*.cpp
#     obj/unittest/support/{driver,random_graph}.o  (minus driver.o for the main binary)
#     obj/subcommand/*.o                # src/subcommand/*.cpp
#     obj/config/allocator_config_X.o  # CONFIG_OBJ
#     -L lib/ <LDFLAGS>
#     lib/libvg.a
#     -lvcflib -lwfa2 ... (LD_LIB_FLAGS, dynamic)
#     [-Wl,-Bstatic]                   # Linux only
#       -lvgio -lhts -ldeflate ...     # LD_STATIC_LIB_FLAGS
#     [-Wl,-Bdynamic]                  # Linux only
#     -lpthread -lm                    # LD_STATIC_LIB_DEPS
#     [libjemalloc.a]                  # LD_EXE_LIB_FLAGS

set(VG_SRC_DIR ${CMAKE_SOURCE_DIR}/src)

# ── Source file globs ──────────────────────────────────────────────────────
file(GLOB SUBCOMMAND_SRCS ${VG_SRC_DIR}/subcommand/*.cpp)
file(GLOB UNITTEST_SRCS   ${VG_SRC_DIR}/unittest/*.cpp)
file(GLOB UNITTEST_SUPPORT_SRCS ${VG_SRC_DIR}/unittest/support/*.cpp)

# driver.cpp provides the main() for per-suite test binaries (Phase 9).
# It is excluded from the main vg binary (uses its own main in main.cpp).
# random_graph.cpp and any other support files DO go in the main binary.
list(REMOVE_ITEM UNITTEST_SUPPORT_SRCS
    ${VG_SRC_DIR}/unittest/support/driver.cpp
)

# ── Allocator config object (conditional) ─────────────────────────────────
# The Makefile conditionally links exactly one config .o based on jemalloc/mimalloc.
if(VG_USE_MIMALLOC)
    add_library(allocator_config OBJECT
        ${VG_SRC_DIR}/config/allocator_config_mimalloc.cpp
    )
    add_dependencies(allocator_config mimalloc_ep)
elseif(VG_USE_JEMALLOC)
    add_library(allocator_config OBJECT
        ${VG_SRC_DIR}/config/allocator_config_jemalloc.cpp
    )
    add_dependencies(allocator_config jemalloc_ep)
else()
    add_library(allocator_config OBJECT
        ${VG_SRC_DIR}/config/allocator_config_system.cpp
    )
endif()

target_include_directories(allocator_config PRIVATE
    ${VG_INC_DIR}
    ${VG_SRC_DIR}
)

# ── vg executable ──────────────────────────────────────────────────────────
add_executable(vg
    ${VG_SRC_DIR}/main.cpp
    ${SUBCOMMAND_SRCS}
    ${UNITTEST_SRCS}
    ${UNITTEST_SUPPORT_SRCS}
    $<TARGET_OBJECTS:allocator_config>
)

set_target_properties(vg PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${VG_BIN_DIR}
)

target_include_directories(vg PRIVATE
    ${VG_INC_DIR}
    ${VG_INC_DIR}/dynamic
    ${CMAKE_SOURCE_DIR}
    ${VG_SRC_DIR}
    ${VG_SRC_DIR}/unittest
    ${VG_SRC_DIR}/unittest/support
    ${VG_SRC_DIR}/subcommand
)

# ── Link order (mirrors Makefile exactly) ─────────────────────────────────
#
# libvg.a → LD_LIB_FLAGS (dynamic) → [Bstatic] LD_STATIC_LIB_FLAGS [Bdynamic]
#         → LD_STATIC_LIB_DEPS → LD_EXE_LIB_FLAGS

# Optional mimalloc pre-link object (PRE_LINK_DEPS)
if(VG_USE_MIMALLOC)
    # Mimalloc must be linked first as a raw object file.
    # Use target_link_options to prepend it to the link command.
    target_link_options(vg PRIVATE
        $<BUILD_INTERFACE:${VG_MIMALLOC_OBJ}>
    )
endif()

# Core link: libvg + dynamic dep libs
target_link_libraries(vg PRIVATE
    # libvg.a (whole static archive — contains all vg object code)
    libvg

    # LD_LIB_FLAGS dynamic section (Makefile order preserved)
    dep_vcflib
    dep_wfa2
    dep_tabixpp
    dep_gssw
    dep_ssw
    dep_sublinearLS
    pthread
    ncurses
    dep_gcsa2
    dep_gbwtgraph
    dep_gbwt
    dep_kff
    dep_divsufsort
    dep_divsufsort64
    dep_libvcfh
    dep_raptor
    dep_pinchesandcacti
    dep_3edgeconnected
    dep_sonlib
    dep_fml
    dep_structures
    dep_libbdsg
    dep_xg
    dep_sdsl
    dep_libhandlegraph
    crypto
)

# Linux-only dynamic link additions (before static section)
if(NOT APPLE)
    target_link_libraries(vg PRIVATE
        dep_dwfl dep_dw dep_dwelf dep_ebl dep_elf
        atomic
    )
endif()

# OpenMP
target_link_libraries(vg PRIVATE ${VG_OMP_LINK_LIBS})

# Boost
target_link_libraries(vg PRIVATE Boost::program_options)

# pkg-config dynamic deps (cairo, libzstd)
target_link_libraries(vg PRIVATE PkgConfig::CAIRO PkgConfig::ZSTD)

# Static section: -Wl,-Bstatic on Linux (empty on Mac)
# These use static linking for better non-PIC code paths (Makefile comment).
target_link_libraries(vg PRIVATE
    $<$<NOT:$<BOOL:${APPLE}>>:-Wl,-Bstatic>

    dep_libvgio
    dep_htslib
    deflate
    z
    bz2
    lzma

    # PKG_CONFIG_STATIC_DEPS: protobuf, jansson (Makefile --static)
    protobuf::libprotobuf
    PkgConfig::JANSSON

    $<$<NOT:$<BOOL:${APPLE}>>:-Wl,-Bdynamic>
)

# LD_STATIC_LIB_DEPS (always dynamic to avoid system-malloc conflicts)
target_link_libraries(vg PRIVATE pthread m)

# LD_EXE_LIB_FLAGS: jemalloc is linked statically (direct path to .a)
if(VG_USE_JEMALLOC)
    target_link_libraries(vg PRIVATE ${VG_LIB_DIR}/libjemalloc.a)
    add_dependencies(vg jemalloc_ep)
endif()

# Library search path: prefer our staged lib/ over system dirs
# (mirrors Makefile's LD_LIB_DIR_FLAGS := -L$(CWD)/lib)
target_link_options(vg PRIVATE
    LINKER:-rpath,${VG_LIB_DIR}
    -L${VG_LIB_DIR}
)
if(NOT APPLE)
    target_link_options(vg PRIVATE LINKER:-rpath,${VG_LIB_DIR})
endif()

# Ensure all dep builds complete before linking vg
add_dependencies(vg libvg)

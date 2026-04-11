# Phase 5: ExternalProject_Add for Makefile/autoconf-based bundled deps
#
# Build order (matches deps/dependency graph from the Makefile):
#   libdeflate → htslib → tabixpp → vcflib
#   sdsl-lite  → gcsa2, gbwt, gbwtgraph (also needs libhandlegraph from Phase 4)
#   sonLib     → pinchesAndCacti
#   ssw, snappy, gssw, libvcfh, fermi-lite, structures, sublinearLS — standalone
#   jemalloc, elfutils — conditional / platform-specific
#   sparsehash — autoconf, headers-only output
#
# At the end of this file we also add:
#   libvgio  (CMake, needs htslib from this file)
#   vcflib   (CMake, needs htslib + tabixpp from this file)
#   libbdsg  (CMake, needs libhandlegraph from Phase 4 + sdsl from this file)

include(ExternalProject)

set(DEPS_DIR ${CMAKE_SOURCE_DIR}/deps)

# Flags propagated into sub-builds
set(VG_SUB_CFLAGS   "-fPIC ${CMAKE_C_FLAGS}")
set(VG_SUB_CXXFLAGS "-fPIC ${CMAKE_CXX_FLAGS}")
# Strip -Xpreprocessor -fopenmp from flags passed to deps that don't like it (snappy, sublinearLS)
string(REPLACE "-Xpreprocessor -fopenmp" "" VG_SUB_CXXFLAGS_NO_OMP "${VG_SUB_CXXFLAGS}")
string(REPLACE "-Werror=return-type"     "" VG_SUB_CXXFLAGS_NO_WERR "${VG_SUB_CXXFLAGS_NO_OMP}")

# On Mac, add AS_INTEGRATED_ASSEMBLER=1 to the make env for asm-heavy deps.
# Use an empty list (not empty string) so ${VG_MAC_ASM_ENV} expands to nothing on Linux.
if(APPLE)
    set(VG_MAC_ASM_ENV AS_INTEGRATED_ASSEMBLER=1)
else()
    set(VG_MAC_ASM_ENV)  # empty list → expands to no arguments
endif()

# ════════════════════════════════════════════════════════════════════════════
# libdeflate  (CMake)
# ════════════════════════════════════════════════════════════════════════════
FetchContent_Declare(
    libdeflate
    GIT_REPOSITORY https://github.com/ebiggers/libdeflate
    GIT_TAG master
)
FetchContent_MakeAvailable(libdeflate)

# ════════════════════════════════════════════════════════════════════════════
# ssw  (standalone Makefile in deps/ssw/src)
# ════════════════════════════════════════════════════════════════════════════
ExternalProject_Add(ssw_ep
    SOURCE_DIR  ${DEPS_DIR}/ssw/src
    BUILD_IN_SOURCE ON
    CONFIGURE_COMMAND ""
    BUILD_COMMAND
        ${CMAKE_MAKE_PROGRAM} clean
        COMMAND ${CMAKE_MAKE_PROGRAM}
            "CFLAGS=${VG_SUB_CFLAGS} -I."
            "CXXFLAGS=${VG_SUB_CXXFLAGS} -I."
    INSTALL_COMMAND
        ar rs ${VG_LIB_DIR}/libssw.a ssw.o ssw_cpp.o
        COMMAND ${CMAKE_COMMAND} -E copy ssw_cpp.h ssw.h ${VG_INC_DIR}
    BUILD_BYPRODUCTS ${VG_LIB_DIR}/libssw.a
)

add_library(dep_ssw STATIC IMPORTED GLOBAL)
set_target_properties(dep_ssw PROPERTIES
    IMPORTED_LOCATION             ${VG_LIB_DIR}/libssw.a
    INTERFACE_INCLUDE_DIRECTORIES ${VG_INC_DIR}
)
add_dependencies(dep_ssw ssw_ep)

# ════════════════════════════════════════════════════════════════════════════
# sdsl-lite  (custom install.sh script; produces libsdsl, libdivsufsort)
# ════════════════════════════════════════════════════════════════════════════
set(CMAKE_POLICY_VERSION_MINIMUM 3.5)
add_subdirectory(${DEPS_DIR}/sdsl-lite ${CMAKE_BINARY_DIR}/build/sdsl-lite EXCLUDE_FROM_ALL)

# Stage sdsl headers into the shared include dir so Makefile-based deps can find them
add_custom_target(sdsl_stage_headers
    COMMAND ${CMAKE_COMMAND} -E copy
        ${CMAKE_BINARY_DIR}/build/sdsl-lite/external/libdivsufsort/include/divsufsort.h
        ${VG_INC_DIR}/sdsl/divsufsort.h
    COMMAND ${CMAKE_COMMAND} -E copy
        ${CMAKE_BINARY_DIR}/build/sdsl-lite/external/libdivsufsort/include/divsufsort64.h
        ${VG_INC_DIR}/sdsl/divsufsort64.h
    # Also stage at the top level: sdsl headers use #include "divsufsort.h" (no sdsl/ prefix)
    # and sdsl's own install puts them at <prefix>/include/, not <prefix>/include/sdsl/
    COMMAND ${CMAKE_COMMAND} -E copy
        ${CMAKE_BINARY_DIR}/build/sdsl-lite/external/libdivsufsort/include/divsufsort.h
        ${VG_INC_DIR}/divsufsort.h
    COMMAND ${CMAKE_COMMAND} -E copy
        ${CMAKE_BINARY_DIR}/build/sdsl-lite/external/libdivsufsort/include/divsufsort64.h
        ${VG_INC_DIR}/divsufsort64.h
    COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${DEPS_DIR}/sdsl-lite/include/sdsl
        ${VG_INC_DIR}/sdsl
    COMMENT "Staging sdsl + divsufsort headers -> ${VG_INC_DIR}"
)
add_dependencies(sdsl_stage_headers sdsl)

add_library(dep_divsufsort STATIC IMPORTED GLOBAL)
set_target_properties(dep_divsufsort PROPERTIES
    IMPORTED_LOCATION ${VG_LIB_DIR}/libdivsufsort.a
    INTERFACE_INCLUDE_DIRECTORIES ${VG_INC_DIR}
)
add_dependencies(sdsl_stage_headers divsufsort)

add_library(dep_divsufsort64 STATIC IMPORTED GLOBAL)
set_target_properties(dep_divsufsort64 PROPERTIES
    IMPORTED_LOCATION ${VG_LIB_DIR}/libdivsufsort64.a
    INTERFACE_INCLUDE_DIRECTORIES ${VG_INC_DIR}
)
add_dependencies(sdsl_stage_headers divsufsort64)

add_library(dep_sdsl STATIC IMPORTED GLOBAL)
set_target_properties(dep_sdsl PROPERTIES
    IMPORTED_LOCATION             ${VG_LIB_DIR}/libsdsl.a
    INTERFACE_INCLUDE_DIRECTORIES ${VG_INC_DIR}
)
add_dependencies(dep_sdsl sdsl_stage_headers)

# ════════════════════════════════════════════════════════════════════════════
# htslib  (autoconf; depends on libdeflate)
# ════════════════════════════════════════════════════════════════════════════
# Note: Makefile runs autoreconf, then guesses --host via a configure dry run.
# For a native build, omitting --host lets autoconf detect it automatically.
ExternalProject_Add(htslib_ep
    SOURCE_DIR  ${DEPS_DIR}/htslib
    BUILD_IN_SOURCE ON
    CONFIGURE_COMMAND
        autoreconf -i
        COMMAND autoheader
        COMMAND autoconf
        COMMAND ./configure
            "CFLAGS=-fPIC -I${DEPS_DIR}/htslib -isystem ${DEPS_DIR}/htslib -I${VG_INC_DIR} ${CMAKE_C_FLAGS}"
            "LDFLAGS=-L${VG_LIB_DIR} ${CMAKE_EXE_LINKER_FLAGS}"
            --with-libdeflate
            --disable-s3
            --disable-gcs
            --disable-libcurl
            --disable-plugins
            --prefix=${CMAKE_BINARY_DIR}
    BUILD_COMMAND
        ${CMAKE_MAKE_PROGRAM} clean
        COMMAND ${CMAKE_MAKE_PROGRAM}
    INSTALL_COMMAND ${CMAKE_MAKE_PROGRAM} install
    DEPENDS libdeflate::libdeflate_shared
    BUILD_BYPRODUCTS
        ${VG_LIB_DIR}/libhts.a
        ${VG_LIB_DIR}/pkgconfig/htslib.pc
        ${VG_LIB_DIR}/libhts.${VG_SHARED_SUFFIX}
)

add_library(dep_htslib STATIC IMPORTED GLOBAL)
set_target_properties(dep_htslib PROPERTIES
    IMPORTED_LOCATION             ${VG_LIB_DIR}/libhts.a
    INTERFACE_INCLUDE_DIRECTORIES ${VG_INC_DIR}
)
add_dependencies(dep_htslib htslib_ep)

# ════════════════════════════════════════════════════════════════════════════
# snappy  (autoconf; standalone; strips -fopenmp from CXXFLAGS)
# ════════════════════════════════════════════════════════════════════════════
ExternalProject_Add(snappy_ep
    SOURCE_DIR  ${DEPS_DIR}/snappy
    BUILD_IN_SOURCE ON
    CONFIGURE_COMMAND
        ./autogen.sh
        COMMAND ./configure
            "CXXFLAGS=-fPIC ${VG_SUB_CXXFLAGS_NO_OMP}"
            --prefix=${CMAKE_BINARY_DIR}
    BUILD_COMMAND
        ${CMAKE_MAKE_PROGRAM} libsnappy.la
            "CXXFLAGS=-fPIC ${VG_SUB_CXXFLAGS_NO_OMP}"
    INSTALL_COMMAND
        ${CMAKE_COMMAND} -E copy .libs/libsnappy.a ${VG_LIB_DIR}/libsnappy.a
        COMMAND ${CMAKE_COMMAND} -E copy
            snappy-c.h snappy-sinksource.h snappy-stubs-public.h snappy.h
            ${VG_INC_DIR}
    BUILD_BYPRODUCTS ${VG_LIB_DIR}/libsnappy.a
)

add_library(dep_snappy STATIC IMPORTED GLOBAL)
set_target_properties(dep_snappy PROPERTIES
    IMPORTED_LOCATION             ${VG_LIB_DIR}/libsnappy.a
    INTERFACE_INCLUDE_DIRECTORIES ${VG_INC_DIR}
)
add_dependencies(dep_snappy snappy_ep)

# ════════════════════════════════════════════════════════════════════════════
# sparsehash  (autoconf; headers-only output after install)
# ════════════════════════════════════════════════════════════════════════════
ExternalProject_Add(sparsehash_ep
    SOURCE_DIR  ${DEPS_DIR}/sparsehash
    BUILD_IN_SOURCE ON
    CONFIGURE_COMMAND
        ./autogen.sh
        COMMAND ./configure
            "LDFLAGS=-L${VG_LIB_DIR} ${CMAKE_EXE_LINKER_FLAGS}"
            --prefix=${CMAKE_BINARY_DIR}
    BUILD_COMMAND ${CMAKE_MAKE_PROGRAM}
    INSTALL_COMMAND ${CMAKE_MAKE_PROGRAM} install
    BUILD_BYPRODUCTS ${VG_INC_DIR}/sparsehash/sparse_hash_map
)

add_library(dep_sparsehash INTERFACE)
target_include_directories(dep_sparsehash INTERFACE ${VG_INC_DIR})
add_dependencies(dep_sparsehash sparsehash_ep)

# ════════════════════════════════════════════════════════════════════════════
# jemalloc  (autoconf; conditional on VG_USE_JEMALLOC; standalone)
# ════════════════════════════════════════════════════════════════════════════
if(VG_USE_JEMALLOC)
    ExternalProject_Add(jemalloc_ep
        SOURCE_DIR  ${DEPS_DIR}/jemalloc
        BUILD_IN_SOURCE ON
        CONFIGURE_COMMAND
            ./autogen.sh
            COMMAND ./configure
                --enable-prof
                --disable-libdl
                --prefix=<SOURCE_DIR>
        BUILD_COMMAND
            ${CMAKE_MAKE_PROGRAM} clean
            COMMAND ${CMAKE_MAKE_PROGRAM}
        INSTALL_COMMAND
            ${CMAKE_COMMAND} -E copy_directory include ${VG_INC_DIR}
            COMMAND ${CMAKE_COMMAND} -E copy lib/libjemalloc.a ${VG_LIB_DIR}/libjemalloc.a
        BUILD_BYPRODUCTS ${VG_LIB_DIR}/libjemalloc.a
    )

    add_library(dep_jemalloc STATIC IMPORTED GLOBAL)
    set_target_properties(dep_jemalloc PROPERTIES
        IMPORTED_LOCATION             ${VG_LIB_DIR}/libjemalloc.a
        INTERFACE_INCLUDE_DIRECTORIES ${VG_INC_DIR}
    )
    add_dependencies(dep_jemalloc jemalloc_ep)
endif()

# ════════════════════════════════════════════════════════════════════════════
# elfutils  (autoconf; Linux-only; multi-step build across sub-directories)
# ════════════════════════════════════════════════════════════════════════════
# The Makefile configures once then builds each sublib individually:
#   libelf, libebl, libdwfl, libdwelf, lib (libeu), libcpu, backends, libdw
if(NOT APPLE)
    ExternalProject_Add(elfutils_ep
        SOURCE_DIR  ${DEPS_DIR}/elfutils
        BUILD_IN_SOURCE ON
        CONFIGURE_COMMAND
            autoreconf -i -f
            COMMAND ./configure
                "CFLAGS=-fPIC ${CMAKE_C_FLAGS}"
                "CXXFLAGS=-fPIC ${CMAKE_CXX_FLAGS}"
                --enable-maintainer-mode
                --disable-libdebuginfod
                --disable-debuginfod
                --prefix=${CMAKE_BINARY_DIR}
        BUILD_COMMAND
            # Build each sub-library individually (as the Makefile does)
            ${CMAKE_MAKE_PROGRAM} -C libelf  clean COMMAND ${CMAKE_MAKE_PROGRAM} -C libelf  libelf.a
            COMMAND ${CMAKE_MAKE_PROGRAM} -C libebl  clean COMMAND ${CMAKE_MAKE_PROGRAM} -C libebl  libebl.a
            COMMAND ${CMAKE_MAKE_PROGRAM} -C libdwfl clean COMMAND ${CMAKE_MAKE_PROGRAM} -C libdwfl libdwfl.a
            COMMAND ${CMAKE_MAKE_PROGRAM} -C libdwelf clean COMMAND ${CMAKE_MAKE_PROGRAM} -C libdwelf libdwelf.a
            COMMAND ${CMAKE_MAKE_PROGRAM} -C lib     clean COMMAND ${CMAKE_MAKE_PROGRAM} -C lib     libeu.a
            COMMAND ${CMAKE_MAKE_PROGRAM} -C libcpu  clean COMMAND ${CMAKE_MAKE_PROGRAM} -C libcpu  libcpu.a
            COMMAND ${CMAKE_MAKE_PROGRAM} -C backends clean "CFLAGS=-fPIC ${CMAKE_C_FLAGS}" COMMAND ${CMAKE_MAKE_PROGRAM} -C backends libebl_backends.a
            COMMAND ${CMAKE_MAKE_PROGRAM} -C libdw   clean "CFLAGS=-fPIC ${CMAKE_C_FLAGS}" COMMAND ${CMAKE_MAKE_PROGRAM} -C libdw   libdw.a known-dwarf.h
        INSTALL_COMMAND
            ${CMAKE_COMMAND} -E make_directory ${VG_INC_DIR}/elfutils
            COMMAND ${CMAKE_COMMAND} -E copy
                libdw/known-dwarf.h libdw/libdw.h libebl/libebl.h
                libelf/elf-knowledge.h version.h libdwfl/libdwfl.h libdwelf/libdwelf.h
                ${VG_INC_DIR}/elfutils
            COMMAND ${CMAKE_COMMAND} -E copy libelf/gelf.h libelf/libelf.h libdw/dwarf.h ${VG_INC_DIR}
            COMMAND ${CMAKE_COMMAND} -E copy
                libebl/libebl.a libdw/libdw.a libdwfl/libdwfl.a
                libdwelf/libdwelf.a libelf/libelf.a
                ${VG_LIB_DIR}
        BUILD_BYPRODUCTS
            ${VG_LIB_DIR}/libelf.a
            ${VG_LIB_DIR}/libebl.a
            ${VG_LIB_DIR}/libdwfl.a
            ${VG_LIB_DIR}/libdwelf.a
            ${VG_LIB_DIR}/libdw.a
    )

    foreach(_elib elf ebl dwfl dwelf dw)
        add_library(dep_${_elib} STATIC IMPORTED GLOBAL)
        set_target_properties(dep_${_elib} PROPERTIES
            IMPORTED_LOCATION             ${VG_LIB_DIR}/lib${_elib}.a
            INTERFACE_INCLUDE_DIRECTORIES ${VG_INC_DIR}
        )
        add_dependencies(dep_${_elib} elfutils_ep)
    endforeach()
endif()

# ════════════════════════════════════════════════════════════════════════════
# gcsa2  (Makefile; depends on sdsl)
# ════════════════════════════════════════════════════════════════════════════
ExternalProject_Add(gcsa2_ep
    SOURCE_DIR  ${DEPS_DIR}/gcsa2
    BUILD_IN_SOURCE ON
    CONFIGURE_COMMAND
        ${CMAKE_COMMAND} -E copy_directory include/gcsa ${VG_INC_DIR}/gcsa
    BUILD_COMMAND
        ${CMAKE_MAKE_PROGRAM} clean
        COMMAND ${CMAKE_MAKE_PROGRAM} directories
        COMMAND ${CMAKE_COMMAND} -E env
            ${VG_MAC_ASM_ENV}
            "CFLAGS=${VG_SUB_CFLAGS}"
            "CXXFLAGS=${VG_SUB_CXXFLAGS}"
            ${CMAKE_MAKE_PROGRAM} lib/libgcsa2.a
                "INC_DIR=${VG_INC_DIR}"
                "LIB_DIR=${VG_LIB_DIR}"
    INSTALL_COMMAND
        ${CMAKE_COMMAND} -E copy lib/libgcsa2.a ${VG_LIB_DIR}/libgcsa2.a
    DEPENDS sdsl_stage_headers
    BUILD_BYPRODUCTS ${VG_LIB_DIR}/libgcsa2.a
)

add_library(dep_gcsa2 STATIC IMPORTED GLOBAL)
set_target_properties(dep_gcsa2 PROPERTIES
    IMPORTED_LOCATION             ${VG_LIB_DIR}/libgcsa2.a
    INTERFACE_INCLUDE_DIRECTORIES ${VG_INC_DIR}
)
add_dependencies(dep_gcsa2 gcsa2_ep)

# ════════════════════════════════════════════════════════════════════════════
# gbwt  (Makefile; depends on sdsl)
# ════════════════════════════════════════════════════════════════════════════
ExternalProject_Add(gbwt_ep
    SOURCE_DIR  ${DEPS_DIR}/gbwt
    BUILD_IN_SOURCE ON
    CONFIGURE_COMMAND
        ${CMAKE_COMMAND} -E copy_directory include/gbwt ${VG_INC_DIR}/gbwt
    BUILD_COMMAND
        ${CMAKE_MAKE_PROGRAM} clean
        COMMAND ${CMAKE_COMMAND} -E env
            ${VG_MAC_ASM_ENV}
            "CFLAGS=${VG_SUB_CFLAGS}"
            "CXXFLAGS=${VG_SUB_CXXFLAGS}"
            ${CMAKE_MAKE_PROGRAM}
                "INC_DIR=${VG_INC_DIR}"
                "LIB_DIR=${VG_LIB_DIR}"
    INSTALL_COMMAND
        ${CMAKE_COMMAND} -E copy lib/libgbwt.a ${VG_LIB_DIR}/libgbwt.a
    DEPENDS sdsl_stage_headers
    BUILD_BYPRODUCTS ${VG_LIB_DIR}/libgbwt.a
)

add_library(dep_gbwt STATIC IMPORTED GLOBAL)
set_target_properties(dep_gbwt PROPERTIES
    IMPORTED_LOCATION             ${VG_LIB_DIR}/libgbwt.a
    INTERFACE_INCLUDE_DIRECTORIES ${VG_INC_DIR}
)
add_dependencies(dep_gbwt gbwt_ep)

# ════════════════════════════════════════════════════════════════════════════
# gbwtgraph  (Makefile; depends on gbwt + sdsl + libhandlegraph)
# ════════════════════════════════════════════════════════════════════════════
ExternalProject_Add(gbwtgraph_ep
    SOURCE_DIR  ${DEPS_DIR}/gbwtgraph
    BUILD_IN_SOURCE ON
    CONFIGURE_COMMAND
        ${CMAKE_COMMAND} -E copy_directory include/gbwtgraph ${VG_INC_DIR}/gbwtgraph
    BUILD_COMMAND
        ${CMAKE_MAKE_PROGRAM} clean
        COMMAND ${CMAKE_COMMAND} -E env
            ${VG_MAC_ASM_ENV}
            "CFLAGS=${VG_SUB_CFLAGS}"
            "CXXFLAGS=${VG_SUB_CXXFLAGS}"
            ${CMAKE_MAKE_PROGRAM}
                "INC_DIR=${VG_INC_DIR}"
                "LIB_DIR=${VG_LIB_DIR}"
    INSTALL_COMMAND
        ${CMAKE_COMMAND} -E copy lib/libgbwtgraph.a ${VG_LIB_DIR}/libgbwtgraph.a
    DEPENDS gbwt_ep sdsl_stage_headers handlegraph_stage_headers
    BUILD_BYPRODUCTS ${VG_LIB_DIR}/libgbwtgraph.a
)

add_library(dep_gbwtgraph STATIC IMPORTED GLOBAL)
set_target_properties(dep_gbwtgraph PROPERTIES
    IMPORTED_LOCATION             ${VG_LIB_DIR}/libgbwtgraph.a
    INTERFACE_INCLUDE_DIRECTORIES ${VG_INC_DIR}
)
add_dependencies(dep_gbwtgraph gbwtgraph_ep)

# ════════════════════════════════════════════════════════════════════════════
# tabixpp  (Makefile; depends on htslib)
# ════════════════════════════════════════════════════════════════════════════
ExternalProject_Add(tabixpp_ep
    SOURCE_DIR  ${DEPS_DIR}/tabixpp
    BUILD_IN_SOURCE ON
    CONFIGURE_COMMAND ""
    BUILD_COMMAND
        ${CMAKE_MAKE_PROGRAM} tabix.o
            "CFLAGS=-fPIC ${CMAKE_C_FLAGS}"
            "CXXFLAGS=-fPIC ${CMAKE_CXX_FLAGS}"
            "INCLUDES=-I${VG_INC_DIR}"
            "HTS_HEADERS="
        COMMAND ar rcs libtabixpp.a tabix.o
    INSTALL_COMMAND
        ${CMAKE_COMMAND} -E copy libtabixpp.a ${VG_LIB_DIR}/libtabixpp.a
        COMMAND ${CMAKE_COMMAND} -E copy tabix.hpp ${VG_INC_DIR}/tabix.hpp
        # Write a minimal pkg-config file for tabixpp (as Makefile does)
        COMMAND ${CMAKE_COMMAND} -P
            ${CMAKE_SOURCE_DIR}/cmake/tabixpp_pkgconfig.cmake
    DEPENDS htslib_ep
    BUILD_BYPRODUCTS ${VG_LIB_DIR}/libtabixpp.a
)

file(WRITE ${CMAKE_SOURCE_DIR}/cmake/tabixpp_pkgconfig.cmake
"file(WRITE \"${VG_PKG_DIR}/tabixpp.pc\"
\"Name: tabixpp
Description: Self-packaged tabixpp
Version: 1.0
Cflags: -I${VG_INC_DIR}
Libs: -L${VG_LIB_DIR} -ltabixpp
\")
")

add_library(dep_tabixpp STATIC IMPORTED GLOBAL)
set_target_properties(dep_tabixpp PROPERTIES
    IMPORTED_LOCATION             ${VG_LIB_DIR}/libtabixpp.a
    INTERFACE_INCLUDE_DIRECTORIES ${VG_INC_DIR}
)
add_dependencies(dep_tabixpp tabixpp_ep)

# ════════════════════════════════════════════════════════════════════════════
# gssw  (Makefile; standalone)
# ════════════════════════════════════════════════════════════════════════════
ExternalProject_Add(gssw_ep
    SOURCE_DIR  ${DEPS_DIR}/gssw
    BUILD_IN_SOURCE ON
    CONFIGURE_COMMAND ""
    BUILD_COMMAND
        ${CMAKE_MAKE_PROGRAM} clean
        COMMAND ${CMAKE_MAKE_PROGRAM}
            "CFLAGS=${VG_SUB_CFLAGS}"
            "CXXFLAGS=${VG_SUB_CXXFLAGS}"
    INSTALL_COMMAND
        ${CMAKE_COMMAND} -E copy lib/libgssw.a ${VG_LIB_DIR}/libgssw.a
        COMMAND ${CMAKE_COMMAND} -E copy src/gssw.h ${VG_INC_DIR}/gssw.h
    BUILD_BYPRODUCTS ${VG_LIB_DIR}/libgssw.a
)

add_library(dep_gssw STATIC IMPORTED GLOBAL)
set_target_properties(dep_gssw PROPERTIES
    IMPORTED_LOCATION             ${VG_LIB_DIR}/libgssw.a
    INTERFACE_INCLUDE_DIRECTORIES ${VG_INC_DIR}
)
add_dependencies(dep_gssw gssw_ep)

# ════════════════════════════════════════════════════════════════════════════
# sonLib  (Makefile; standalone)
# ════════════════════════════════════════════════════════════════════════════
ExternalProject_Add(sonlib_ep
    SOURCE_DIR  ${DEPS_DIR}/sonLib
    BUILD_IN_SOURCE ON
    CONFIGURE_COMMAND ""
    BUILD_COMMAND
        ${CMAKE_MAKE_PROGRAM} clean
        COMMAND ${CMAKE_COMMAND} -E env
            "kyotoTycoonLib="
            "CFLAGS=${VG_SUB_CFLAGS}"
            "CXXFLAGS=${VG_SUB_CXXFLAGS}"
            ${CMAKE_MAKE_PROGRAM}
    INSTALL_COMMAND
        ${CMAKE_COMMAND} -E copy lib/sonLib.a ${VG_LIB_DIR}/libsonlib.a
        COMMAND ${CMAKE_COMMAND} -E make_directory ${VG_INC_DIR}/sonLib
        COMMAND ${CMAKE_COMMAND} -E copy_directory lib ${VG_INC_DIR}/sonLib
        # copy_directory will copy *.h — that's intentional
    BUILD_BYPRODUCTS ${VG_LIB_DIR}/libsonlib.a
)

add_library(dep_sonlib STATIC IMPORTED GLOBAL)
set_target_properties(dep_sonlib PROPERTIES
    IMPORTED_LOCATION             ${VG_LIB_DIR}/libsonlib.a
    INTERFACE_INCLUDE_DIRECTORIES ${VG_INC_DIR}
)
add_dependencies(dep_sonlib sonlib_ep)

# ════════════════════════════════════════════════════════════════════════════
# pinchesAndCacti  (Makefile; depends on sonLib)
# ════════════════════════════════════════════════════════════════════════════
ExternalProject_Add(pinchescacti_ep
    SOURCE_DIR  ${DEPS_DIR}/pinchesAndCacti
    BUILD_IN_SOURCE ON
    CONFIGURE_COMMAND ""
    BUILD_COMMAND
        ${CMAKE_MAKE_PROGRAM} clean
        COMMAND ${CMAKE_COMMAND} -E env
            "CFLAGS=${VG_SUB_CFLAGS}"
            "CXXFLAGS=${VG_SUB_CXXFLAGS}"
            ${CMAKE_MAKE_PROGRAM}
    INSTALL_COMMAND
        # Libraries land in sonLib's lib/ after pinchesAndCacti builds
        ${CMAKE_COMMAND} -E copy
            ${DEPS_DIR}/sonLib/lib/stPinchesAndCacti.a
            ${VG_LIB_DIR}/libpinchesandcacti.a
        COMMAND ${CMAKE_COMMAND} -E copy
            ${DEPS_DIR}/sonLib/lib/3EdgeConnected.a
            ${VG_LIB_DIR}/lib3edgeconnected.a
        COMMAND ${CMAKE_COMMAND} -E make_directory ${VG_INC_DIR}/sonLib
        COMMAND ${CMAKE_COMMAND} -E copy_directory
            ${DEPS_DIR}/sonLib/lib ${VG_INC_DIR}/sonLib
    DEPENDS sonlib_ep
    BUILD_BYPRODUCTS
        ${VG_LIB_DIR}/libpinchesandcacti.a
        ${VG_LIB_DIR}/lib3edgeconnected.a
)

add_library(dep_pinchesandcacti STATIC IMPORTED GLOBAL)
set_target_properties(dep_pinchesandcacti PROPERTIES
    IMPORTED_LOCATION             ${VG_LIB_DIR}/libpinchesandcacti.a
)
add_dependencies(dep_pinchesandcacti pinchescacti_ep)

add_library(dep_3edgeconnected STATIC IMPORTED GLOBAL)
set_target_properties(dep_3edgeconnected PROPERTIES
    IMPORTED_LOCATION             ${VG_LIB_DIR}/lib3edgeconnected.a
)
add_dependencies(dep_3edgeconnected pinchescacti_ep)

# ════════════════════════════════════════════════════════════════════════════
# libVCFH  (Makefile; standalone)
# ════════════════════════════════════════════════════════════════════════════
ExternalProject_Add(libvcfh_ep
    SOURCE_DIR  ${DEPS_DIR}/libVCFH
    BUILD_IN_SOURCE ON
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ${CMAKE_MAKE_PROGRAM}
    INSTALL_COMMAND
        ${CMAKE_COMMAND} -E copy libvcfh.a    ${VG_LIB_DIR}/libvcfh.a
        COMMAND ${CMAKE_COMMAND} -E copy vcfheader.hpp ${VG_INC_DIR}/vcfheader.hpp
    BUILD_BYPRODUCTS ${VG_LIB_DIR}/libvcfh.a
)

add_library(dep_libvcfh STATIC IMPORTED GLOBAL)
set_target_properties(dep_libvcfh PROPERTIES
    IMPORTED_LOCATION             ${VG_LIB_DIR}/libvcfh.a
    INTERFACE_INCLUDE_DIRECTORIES ${VG_INC_DIR}
)
add_dependencies(dep_libvcfh libvcfh_ep)

# ════════════════════════════════════════════════════════════════════════════
# fermi-lite  (Makefile; standalone)
# ════════════════════════════════════════════════════════════════════════════
ExternalProject_Add(fml_ep
    SOURCE_DIR  ${DEPS_DIR}/fermi-lite
    BUILD_IN_SOURCE ON
    CONFIGURE_COMMAND ""
    BUILD_COMMAND
        ${CMAKE_MAKE_PROGRAM} clean
        COMMAND ${CMAKE_MAKE_PROGRAM}
            "CFLAGS=${VG_SUB_CFLAGS}"
            "CXXFLAGS=${VG_SUB_CXXFLAGS}"
    INSTALL_COMMAND
        ${CMAKE_COMMAND} -E copy libfml.a ${VG_LIB_DIR}/libfml.a
        COMMAND ${CMAKE_COMMAND} -E copy_directory . ${VG_INC_DIR}
        # copies all *.h from fermi-lite to include/
    BUILD_BYPRODUCTS ${VG_LIB_DIR}/libfml.a
)

add_library(dep_fml STATIC IMPORTED GLOBAL)
set_target_properties(dep_fml PROPERTIES
    IMPORTED_LOCATION             ${VG_LIB_DIR}/libfml.a
    INTERFACE_INCLUDE_DIRECTORIES ${VG_INC_DIR}
)
add_dependencies(dep_fml fml_ep)

# ════════════════════════════════════════════════════════════════════════════
# structures  (Makefile; standalone)
# ════════════════════════════════════════════════════════════════════════════
ExternalProject_Add(structures_ep
    SOURCE_DIR  ${DEPS_DIR}/structures
    BUILD_IN_SOURCE ON
    CONFIGURE_COMMAND ""
    BUILD_COMMAND
        ${CMAKE_MAKE_PROGRAM} clean
        COMMAND ${CMAKE_MAKE_PROGRAM} lib/libstructures.a
    INSTALL_COMMAND
        ${CMAKE_COMMAND} -E copy lib/libstructures.a ${VG_LIB_DIR}/libstructures.a
        COMMAND ${CMAKE_COMMAND} -E copy_directory src/include/structures ${VG_INC_DIR}/structures
    BUILD_BYPRODUCTS ${VG_LIB_DIR}/libstructures.a
)

add_library(dep_structures STATIC IMPORTED GLOBAL)
set_target_properties(dep_structures PROPERTIES
    IMPORTED_LOCATION             ${VG_LIB_DIR}/libstructures.a
    INTERFACE_INCLUDE_DIRECTORIES ${VG_INC_DIR}
)
add_dependencies(dep_structures structures_ep)

# ════════════════════════════════════════════════════════════════════════════
# sublinear-Li-Stephens  (Makefile; uses htslib via INCLUDE_FLAGS)
# ════════════════════════════════════════════════════════════════════════════
# Must strip -Werror=return-type (upstream bug workaround per Makefile comment)
ExternalProject_Add(sublinearls_ep
    SOURCE_DIR  ${DEPS_DIR}/sublinear-Li-Stephens
    BUILD_IN_SOURCE ON
    CONFIGURE_COMMAND ""
    BUILD_COMMAND
        ${CMAKE_MAKE_PROGRAM} clean
        COMMAND ${CMAKE_MAKE_PROGRAM} libs
            "CFLAGS=-fPIC ${VG_SUB_CXXFLAGS_NO_WERR}"
            "CXXFLAGS=-fPIC ${VG_SUB_CXXFLAGS_NO_WERR}"
            "INCLUDE_FLAGS=-I${VG_INC_DIR}"
    INSTALL_COMMAND
        ${CMAKE_COMMAND} -E copy lib/libsublinearLS.a ${VG_LIB_DIR}/libsublinearLS.a
        COMMAND ${CMAKE_COMMAND} -E make_directory ${VG_INC_DIR}/sublinearLS
        COMMAND ${CMAKE_COMMAND} -E copy_directory src ${VG_INC_DIR}/sublinearLS
    DEPENDS htslib_ep
    BUILD_BYPRODUCTS ${VG_LIB_DIR}/libsublinearLS.a
)

add_library(dep_sublinearLS STATIC IMPORTED GLOBAL)
set_target_properties(dep_sublinearLS PROPERTIES
    IMPORTED_LOCATION             ${VG_LIB_DIR}/libsublinearLS.a
    INTERFACE_INCLUDE_DIRECTORIES ${VG_INC_DIR}
)
add_dependencies(dep_sublinearLS sublinearls_ep)

# ════════════════════════════════════════════════════════════════════════════
# vcflib  (CMake; depends on htslib + tabixpp — placed here because of deps)
# ════════════════════════════════════════════════════════════════════════════
# Makefile: cmake ... --target vcflib vcf2tsv wfa2_static
ExternalProject_Add(vcflib_ep
    SOURCE_DIR  ${DEPS_DIR}/vcflib
    BINARY_DIR  ${CMAKE_BINARY_DIR}/build/vcflib
    CMAKE_ARGS
        -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
        -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
        -DCMAKE_VERBOSE_MAKEFILE=ON
        -DZIG=OFF
        -DCMAKE_C_FLAGS=${CMAKE_C_FLAGS}
        -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}
        -DCMAKE_BUILD_TYPE=RelWithDebInfo
        -DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}
        # CMAKE_PREFIX_PATH: ${CMAKE_BINARY_DIR} makes pkg_check_modules find our
        # installed htslib.pc / tabixpp.pc in ${CMAKE_BINARY_DIR}/lib/pkgconfig
        "-DCMAKE_PREFIX_PATH=${CMAKE_BINARY_DIR}\;/usr\;${VG_OMP_PREFIXES_STR}"
        "-DPYTHON_EXECUTABLE=${Python3_EXECUTABLE}"
    BUILD_COMMAND
        ${CMAKE_MAKE_PROGRAM} vcflib vcf2tsv wfa2_static
    INSTALL_COMMAND
        ${CMAKE_COMMAND} -E copy_directory
            ${DEPS_DIR}/vcflib/contrib/filevercmp    ${VG_INC_DIR}
        COMMAND ${CMAKE_COMMAND} -E copy_directory
            ${DEPS_DIR}/vcflib/contrib/fastahack     ${VG_INC_DIR}
        COMMAND ${CMAKE_COMMAND} -E copy_directory
            ${DEPS_DIR}/vcflib/contrib/smithwaterman ${VG_INC_DIR}
        COMMAND ${CMAKE_COMMAND} -E copy
            ${DEPS_DIR}/vcflib/contrib/intervaltree/IntervalTree.h ${VG_INC_DIR}
        COMMAND ${CMAKE_COMMAND} -E copy_directory
            ${DEPS_DIR}/vcflib/contrib/multichoose   ${VG_INC_DIR}
        COMMAND ${CMAKE_COMMAND} -E copy_directory
            ${DEPS_DIR}/vcflib/src                   ${VG_INC_DIR}
        COMMAND ${CMAKE_COMMAND} -E copy
            <BINARY_DIR>/libvcflib.a                 ${VG_LIB_DIR}/libvcflib.a
        COMMAND ${CMAKE_COMMAND} -E copy
            <BINARY_DIR>/contrib/WFA2-lib/libwfa2.a  ${VG_LIB_DIR}/libwfa2.a
    DEPENDS htslib_ep tabixpp_ep
    BUILD_BYPRODUCTS
        ${VG_LIB_DIR}/libvcflib.a
        ${VG_LIB_DIR}/libwfa2.a
)

add_library(dep_vcflib STATIC IMPORTED GLOBAL)
set_target_properties(dep_vcflib PROPERTIES
    IMPORTED_LOCATION             ${VG_LIB_DIR}/libvcflib.a
    INTERFACE_INCLUDE_DIRECTORIES ${VG_INC_DIR}
)
add_dependencies(dep_vcflib vcflib_ep)

add_library(dep_wfa2 STATIC IMPORTED GLOBAL)
set_target_properties(dep_wfa2 PROPERTIES
    IMPORTED_LOCATION             ${VG_LIB_DIR}/libwfa2.a
)
add_dependencies(dep_wfa2 vcflib_ep)

# ════════════════════════════════════════════════════════════════════════════
# libvgio  (CMake; depends on htslib + libhandlegraph + protobuf)
# ════════════════════════════════════════════════════════════════════════════
add_subdirectory(${DEPS_DIR}/libvgio ${CMAKE_BINARY_DIR}/build/libvgio)
add_custom_target(vgio_stage_headers
    COMMAND 
        mkdir -p ${VG_INC_DIR}/vg/io
    COMMAND 
        protoc -I${DEPS_DIR}/libvgio/deps/
               --cpp_out=${VG_INC_DIR}/vg
               ${DEPS_DIR}/libvgio/deps/vg.proto
    COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${DEPS_DIR}/libvgio/include/vg/io
        ${VG_INC_DIR}/vg/io
)
add_dependencies(vgio_stage_headers vgio)

add_library(dep_libvgio STATIC IMPORTED GLOBAL)
set_target_properties(dep_libvgio PROPERTIES
    IMPORTED_LOCATION             ${VG_LIB_DIR}/libvgio.a
    INTERFACE_INCLUDE_DIRECTORIES ${VG_INC_DIR}
)
add_dependencies(dep_libvgio vgio_stage_headers)
target_link_libraries(dep_libvgio INTERFACE dep_htslib dep_libhandlegraph)

# ════════════════════════════════════════════════════════════════════════════
# libbdsg  (CMake; depends on libhandlegraph + sdsl + sparsepp + dynamic + mio)
# ════════════════════════════════════════════════════════════════════════════
set(RUN_DOXYGEN OFF CACHE BOOL "Whether to run Doxygen for libbdsg")
add_subdirectory(${DEPS_DIR}/libbdsg ${CMAKE_BINARY_DIR}/build/libbdsg)
add_custom_target(bdsg_stage_headers
    COMMAND   
        mkdir -p ${VG_INC_DIR}/bdsg
    COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${DEPS_DIR}/libbdsg/bdsg/include/bdsg
        ${VG_INC_DIR}/bdsg
)
# libbdsg's CMakeLists only builds a SHARED lib; create a static archive from its
# OBJECT library (bdsg_objs) so dep_libbdsg's IMPORTED_LOCATION is satisfied.
add_library(libbdsg_archive STATIC $<TARGET_OBJECTS:bdsg_objs>)
set_target_properties(libbdsg_archive PROPERTIES
    OUTPUT_NAME               bdsg
    ARCHIVE_OUTPUT_DIRECTORY  ${VG_LIB_DIR}
    POSITION_INDEPENDENT_CODE ON
)
add_library(dep_libbdsg STATIC IMPORTED GLOBAL)
set_target_properties(dep_libbdsg PROPERTIES
    IMPORTED_LOCATION             ${VG_LIB_DIR}/libbdsg.a
    INTERFACE_INCLUDE_DIRECTORIES ${VG_INC_DIR}
)
target_link_libraries(dep_libbdsg INTERFACE libbdsg_archive dep_libhandlegraph dep_sparsehash)
add_dependencies(dep_libbdsg bdsg_stage_headers libbdsg_archive)
cmake_minimum_required(VERSION 3.22)

project(vg
    LANGUAGES CXX C
    DESCRIPTION "Variation Graph Tools"
)

# ---------------------------------------------------------------------------
# Project setup
# ---------------------------------------------------------------------------
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POLICY_VERSION_MINIMUM 3.5)

set(CMAKE_INSTALL_INCLUDEDIR include CACHE PATH "" FORCE)
set(CMAKE_INSTALL_LIBDIR     lib     CACHE PATH "" FORCE)
set(BUILD_PYTHON_BINDINGS OFF CACHE PATH "" FORCE)
set(RUN_DOXYGEN OFF CACHE PATH "" FORCE)

include(FetchContent)
include(ExternalProject)
include(GNUInstallDirs)

# Build options. Allocator choices are mutually exclusive; jemalloc wins by default.
option(VG_USE_JEMALLOC      "Use jemalloc as the allocator (default)"       ON)
option(VG_USE_MIMALLOC      "Use mimalloc as the allocator"                 OFF)
option(VG_USE_SYSTEM_ALLOC  "Use the system libc allocator"                 OFF)
option(VG_BUILD_UNITTESTS   "Compile and link the in-tree Catch2 tests"     ON)
option(VG_USE_BACKWARD      "Link backward-cpp for stack traces on crash"   ON)

if(VG_USE_MIMALLOC OR VG_USE_SYSTEM_ALLOC)
    set(VG_USE_JEMALLOC OFF)
endif()

set(VG_DEPS_DIR  "${CMAKE_SOURCE_DIR}/deps")
set(VG_STAGE_DIR "${CMAKE_BINARY_DIR}/external")

# OpenMP is non-negotiable for vg.
find_package(OpenMP REQUIRED)
find_package(Threads REQUIRED)
find_package(ZLIB    REQUIRED)
find_package(BZip2   REQUIRED)
find_package(LibLZMA REQUIRED)
find_package(CURL)

# ---------------------------------------------------------------------------
# 2a. FetchContent: deps that ship a CMakeLists.txt and have a canonical
#     upstream we already trust. The on-disk deps/<name>/ checkouts of these
#     are intentionally ignored — they are untracked working-tree leftovers.
# ---------------------------------------------------------------------------
FetchContent_Declare(sdsl-lite
    GIT_REPOSITORY https://github.com/vgteam/sdsl-lite.git
    GIT_TAG        5c425837929b6f323fd7cccfa6c15ffcec6d94bb
    SOURCE_DIR     ${VG_DEPS_DIR}/sdsl-lite)
FetchContent_Declare(libhandlegraph
    GIT_REPOSITORY https://github.com/vgteam/libhandlegraph.git
    GIT_TAG        ae42926da7743a7108db0474bfcf1a1cdac34abb
    SOURCE_DIR     ${VG_DEPS_DIR}/libhandlegraph)
FetchContent_MakeAvailable(
    libhandlegraph
    sdsl-lite
)
target_include_directories(sdsl INTERFACE
    $<BUILD_INTERFACE:${sdsl-lite_BINARY_DIR}/external/libdivsufsort/include>
    $<BUILD_INTERFACE:${sdsl-lite_BINARY_DIR}/external/libdivsufsort/include/64>)

FetchContent_Declare(libvgio
    GIT_REPOSITORY https://github.com/vgteam/libvgio.git
    GIT_TAG        c23b62300c2f9274fc7302a7febc45bf3a78890a
    EXCLUDE_FROM_ALL
)
FetchContent_Declare(libbdsg
    GIT_REPOSITORY https://github.com/vgteam/libbdsg.git
    GIT_TAG        789ad6daf74192b477fa49fb066a8cf0df76a16d
    EXCLUDE_FROM_ALL
)
FetchContent_Declare(Complete-Striped-Smith-Waterman-Library
    GIT_REPOSITORY https://github.com/vgteam/Complete-Striped-Smith-Waterman-Library.git
    GIT_TAG        f53d459841b632ac6ab0e59e1a0ef6d056034631
    EXCLUDE_FROM_ALL
)
FetchContent_Declare(gssw
    GIT_REPOSITORY https://github.com/vgteam/gssw.git
    GIT_TAG        594357fd2473342e6537979a091c1c6d1480b953
    EXCLUDE_FROM_ALL
)
FetchContent_Declare(xg
    GIT_REPOSITORY https://github.com/vgteam/xg.git
    GIT_TAG        vg-cmake-hublabel
    EXCLUDE_FROM_ALL
)

# uses Makefile
FetchContent_Declare(cpp_progress_bar
    GIT_REPOSITORY https://github.com/vgteam/cpp_progress_bar.git
    GIT_TAG        487bb7c9c5f028daaf3b4e24b51d9e783517768f)
# uses Makefile
FetchContent_Declare(dozeu_
    GIT_REPOSITORY https://github.com/vgteam/dozeu.git
    GIT_TAG        999a0ed5457e5a2ab3f47bb62b00fd3355b66dbf
    SOURCE_DIR     ${CMAKE_BINARY_DIR}/deps/dozeu_include/dozeu
)

set(CMAKE_DISABLE_FIND_PACKAGE_pybind11 ON CACHE BOOL "" FORCE)
set(ZIG OFF CACHE BOOL "" FORCE)
FetchContent_Declare(vcflib
    GIT_REPOSITORY https://github.com/vgteam/vcflib.git
    GIT_TAG        pybind-optional
    EXCLUDE_FROM_ALL
)
FetchContent_Declare(structures
    GIT_REPOSITORY https://github.com/vgteam/structures.git
    GIT_TAG        0b31cacac76d55d30f6f653950e8c74aeb8896f6)


FetchContent_MakeAvailable(
    libvgio
    libbdsg
    Complete-Striped-Smith-Waterman-Library
    gssw
    cpp_progress_bar
    dozeu_
    vcflib
    structures
)

FetchContent_GetProperties(xg)
if(NOT xg_POPULATED)
    FetchContent_Populate(xg)   # downloads sources, does NOT add_subdirectory
endif()

# Patch upstream gssw which doesn't expose its include dir.
target_include_directories(libgssw INTERFACE ${gssw_SOURCE_DIR}/src)

# Patch upstream vcflib which doesn't expose its include dir.
target_include_directories(vcflib INTERFACE 
    ${vcflib_SOURCE_DIR}/src 
    ${vcflib_SOURCE_DIR}/contrib/filevercmp 
    ${vcflib_SOURCE_DIR}/contrib/fastahack
    ${vcflib_SOURCE_DIR}/contrib/intervaltree
    ${vcflib_SOURCE_DIR}/contrib/smithwaterman
    ${vcflib_SOURCE_DIR}/contrib/multichoose
)

# external CMake
FetchContent_Declare(backward-cpp
    GIT_REPOSITORY https://github.com/bombela/backward-cpp.git
    GIT_TAG        0bfd0a07a61551413ccd2ab9a9099af3bad40681)
FetchContent_Declare(raptor
    GIT_REPOSITORY https://github.com/dajobe/raptor.git
    GIT_TAG        db8e1a022f2198f207dafe1110344d045d0f3e06
    EXCLUDE_FROM_ALL
)

FetchContent_MakeAvailable(
    backward-cpp
    raptor
)
 # Expose raptor2 headers for consumers; BUILD_INTERFACE only, no install needed.
 target_include_directories(raptor2 INTERFACE
     $<BUILD_INTERFACE:${raptor_BINARY_DIR}/src>
     $<BUILD_INTERFACE:${raptor_SOURCE_DIR}/src>
 )

 # raptor's utils/CMakeLists.txt uses ${CMAKE_BINARY_DIR}/src for raptor_config.h,
 # which is wrong when raptor is a FetchContent dependency. Its generated headers
 # live in ${raptor_BINARY_DIR}/src instead.
 foreach(_t rapper rdfdiff)
     if(TARGET ${_t})
         target_include_directories(${_t} PRIVATE "${raptor_BINARY_DIR}/src")
     endif()
 endforeach()
target_compile_options(turtle_lexer_test PRIVATE -Wno-incompatible-pointer-types)

# ---------------------------------------------------------------------------
# 2b. Vendored deps: small CMake-native libs that exist only locally
# ---------------------------------------------------------------------------
# These all have a CMakeLists.txt under deps/<name>/ — pull them in directly.
foreach(dep IN ITEMS atomic_queue DYNAMIC mio kff-cpp-api mmmultimap)
    if(EXISTS "${VG_DEPS_DIR}/${dep}/CMakeLists.txt")
        add_subdirectory("${VG_DEPS_DIR}/${dep}" "${CMAKE_BINARY_DIR}/deps/${dep}" EXCLUDE_FROM_ALL)
    endif()
endforeach()

# ---------------------------------------------------------------------------
# 2b'. Vendored sources wrapped as add_library (no upstream build system, or
#      upstream is just a Makefile we don't want to invoke).
# ---------------------------------------------------------------------------

# ssw — single C TU + a C++ wrapper header
add_library(vg_ssw STATIC
    ${VG_DEPS_DIR}/ssw/src/ssw.c)
target_include_directories(vg_ssw PUBLIC ${VG_DEPS_DIR}/ssw/src)

# sha1
add_library(vg_sha1 STATIC ${VG_DEPS_DIR}/sha1/sha1.cpp)
target_include_directories(vg_sha1 PUBLIC ${VG_DEPS_DIR}/sha1)

# fermi-lite — C-only de Bruijn assembler. Glob its top-level .c files,
# excluding the example/main programs.
file(GLOB FERMI_LITE_SOURCES CONFIGURE_DEPENDS ${VG_DEPS_DIR}/fermi-lite/*.c)
list(FILTER FERMI_LITE_SOURCES EXCLUDE REGEX "(example|fml-asm)\\.c$")
add_library(vg_fermi_lite STATIC ${FERMI_LITE_SOURCES})
target_include_directories(vg_fermi_lite PUBLIC ${VG_DEPS_DIR}/fermi-lite)
target_link_libraries(vg_fermi_lite PUBLIC ZLIB::ZLIB)

# libdeflate — small C library; wrap its lib/ directory rather than invoking
# its Makefile. Header lives at top of deps/libdeflate.
file(GLOB_RECURSE LIBDEFLATE_SOURCES CONFIGURE_DEPENDS
    ${VG_DEPS_DIR}/libdeflate/lib/*.c)
add_library(vg_libdeflate STATIC ${LIBDEFLATE_SOURCES})
target_include_directories(vg_libdeflate
    PUBLIC  ${VG_DEPS_DIR}/libdeflate
    PRIVATE ${VG_DEPS_DIR}/libdeflate/common
            ${VG_DEPS_DIR}/libdeflate/lib)
target_compile_definitions(vg_libdeflate PRIVATE LIBDEFLATE_STATIC)

# xg — compile xg.cpp directly. Don't use upstream xg's CMakeLists (it drags in
 # conflicting FetchContent for handlegraph/sdsl and broken macOS OpenMP).
 add_library(vg_xg STATIC ${xg_SOURCE_DIR}/src/xg.cpp)
 target_include_directories(vg_xg PUBLIC
     ${xg_SOURCE_DIR}/src
     ${xg_SOURCE_DIR}/deps/ips4o
     ${xg_SOURCE_DIR}/deps/atomic_queue/include/atomic_queue   
     ${xg_SOURCE_DIR}/deps/mmmulti/src
     ${xg_SOURCE_DIR}/deps/gfakluge/src
     ${xg_SOURCE_DIR}/deps/gfakluge/src/tinyFA)
 target_include_directories(vg_xg PRIVATE
     ${sdsl-lite_SOURCE_DIR}/include
     ${sdsl-lite_BINARY_DIR}/external/libdivsufsort/include
     ${sdsl-lite_BINARY_DIR}/external/libdivsufsort/include/64
     ${libhandlegraph_SOURCE_DIR}/src/include)
 target_link_libraries(vg_xg PUBLIC
     mio handlegraph_static sdsl divsufsort divsufsort64
     OpenMP::OpenMP_CXX Threads::Threads ZLIB::ZLIB)
 target_compile_definitions(vg_xg PRIVATE
     _LIBCPP_ENABLE_CXX17_REMOVED_RESULT_OF=1
     _LIBCPP_ENABLE_CXX20_REMOVED_TYPE_TRAITS=1)
 target_compile_options(vg_xg PRIVATE -Wno-missing-template-arg-list-after-template-kw)

# tabixpp — thin C++ wrapper over htslib (htslib comes from ExternalProject)
add_library(vg_tabixpp STATIC ${VG_DEPS_DIR}/tabixpp/tabix.cpp)
target_include_directories(vg_tabixpp PUBLIC ${VG_DEPS_DIR}/tabixpp)

# libVCFH
add_library(vg_libvcfh STATIC ${VG_DEPS_DIR}/libVCFH/vcfheader.cpp)
target_include_directories(vg_libvcfh PUBLIC ${VG_DEPS_DIR}/libVCFH)

# sublinear-Li-Stephens — has many .cpp files in src/, no upstream CMake we
# want. Compile as static lib.
file(GLOB SUBLINEAR_SOURCES CONFIGURE_DEPENDS
    ${VG_DEPS_DIR}/sublinear-Li-Stephens/src/*.cpp)
list(FILTER SUBLINEAR_SOURCES EXCLUDE REGEX "(linhapexample)\\.")
if(SUBLINEAR_SOURCES)
    add_library(vg_sublinear STATIC ${SUBLINEAR_SOURCES})
    target_include_directories(vg_sublinear PUBLIC ${VG_DEPS_DIR}/sublinear-Li-Stephens/src)
endif()

# cpp_progress_bar
add_library(progress_bar STATIC ${cpp_progress_bar_SOURCE_DIR}/progress_bar.cpp)
target_include_directories(progress_bar PUBLIC ${cpp_progress_bar_SOURCE_DIR}/)

# dozeu
add_library(dozeu INTERFACE)
target_include_directories(dozeu INTERFACE ${CMAKE_BINARY_DIR}/deps/dozeu_include/)

# Header-only INTERFACE libraries
add_library(vg_headers_only INTERFACE)
target_include_directories(vg_headers_only INTERFACE
    ${VG_DEPS_DIR}/ips4o
    ${VG_DEPS_DIR}/lru_cache
    ${VG_DEPS_DIR}/BBHash
    ${VG_DEPS_DIR}/sparsepp
    ${VG_DEPS_DIR}/sparsehash/src      # uses local config.h; falls back to system if missing
    ${VG_DEPS_DIR}/sparsehash/src/google
)

# ---------------------------------------------------------------------------
# 2c. ExternalProject: deps with autotools or custom Makefiles
# ---------------------------------------------------------------------------
# htslib (autotools)
set(HTSLIB_INSTALL_DIR "${VG_STAGE_DIR}/htslib")
ExternalProject_Add(htslib_ep
    SOURCE_DIR        ${VG_DEPS_DIR}/htslib/htslib
    BUILD_IN_SOURCE   1
    CONFIGURE_COMMAND autoreconf -i COMMAND ./configure --prefix=${HTSLIB_INSTALL_DIR} --disable-bz2 --disable-lzma
    BUILD_COMMAND     $(MAKE)
    INSTALL_COMMAND   $(MAKE) install
    BUILD_BYPRODUCTS  ${HTSLIB_INSTALL_DIR}/lib/libhts.a)
file(MAKE_DIRECTORY ${HTSLIB_INSTALL_DIR}/include)
add_library(vg_htslib STATIC IMPORTED)
set_target_properties(vg_htslib PROPERTIES
    IMPORTED_LOCATION ${HTSLIB_INSTALL_DIR}/lib/libhts.a
    INTERFACE_INCLUDE_DIRECTORIES ${HTSLIB_INSTALL_DIR}/include)
add_dependencies(vg_htslib htslib_ep)
target_link_libraries(vg_htslib INTERFACE ZLIB::ZLIB Threads::Threads)

# Make tabixpp's #include <htslib/...> resolve.
target_link_libraries(vg_tabixpp PUBLIC vg_htslib)
# ditto for sublinear
target_link_libraries(vg_sublinear PUBLIC vg_htslib)

# ---------------------------------------------------------------------------
# Unified staging prefix for the Makefile-based deps (gbwt/gcsa2/gbwtgraph).
# Their hand-written Makefiles + sdsl-lite's Make.helper expect a single
# include/ + lib/ tree containing sdsl, divsufsort{,64}, handlegraph, and
# (for gbwtgraph) gbwt headers. We build that tree here.
# ---------------------------------------------------------------------------
set(MAKEDEPS_PREFIX     "${VG_STAGE_DIR}/makedeps")
set(MAKEDEPS_INCLUDE    "${MAKEDEPS_PREFIX}/include")
set(MAKEDEPS_LIB        "${MAKEDEPS_PREFIX}/lib")
file(MAKE_DIRECTORY ${MAKEDEPS_INCLUDE} ${MAKEDEPS_LIB})

# Pre-stage: sdsl headers, divsufsort headers (from sdsl-lite's build tree),
# handlegraph headers, and the static libs that gbwt/gcsa2 link against.
add_custom_command(
    OUTPUT  ${MAKEDEPS_INCLUDE}/.stamp
    COMMAND ${CMAKE_COMMAND} -E make_directory ${MAKEDEPS_INCLUDE} ${MAKEDEPS_LIB}
    COMMAND ${CMAKE_COMMAND} -E copy_directory
            ${sdsl-lite_SOURCE_DIR}/include ${MAKEDEPS_INCLUDE}
    COMMAND ${CMAKE_COMMAND} -E copy_directory
            ${sdsl-lite_BINARY_DIR}/external/libdivsufsort/include ${MAKEDEPS_INCLUDE}
    COMMAND ${CMAKE_COMMAND} -E copy_directory
            ${libhandlegraph_SOURCE_DIR}/src/include ${MAKEDEPS_INCLUDE}
    COMMAND ${CMAKE_COMMAND} -E copy
            $<TARGET_FILE:sdsl> ${MAKEDEPS_LIB}/libsdsl.a
    COMMAND ${CMAKE_COMMAND} -E copy
            $<TARGET_FILE:divsufsort>   ${MAKEDEPS_LIB}/libdivsufsort.a
    COMMAND ${CMAKE_COMMAND} -E copy
            $<TARGET_FILE:divsufsort64> ${MAKEDEPS_LIB}/libdivsufsort64.a
    COMMAND ${CMAKE_COMMAND} -E touch ${MAKEDEPS_INCLUDE}/.stamp
    DEPENDS sdsl divsufsort divsufsort64 handlegraph_static
    COMMENT "Staging sdsl + divsufsort + handlegraph for Makefile-based deps")
add_custom_target(makedeps_stage DEPENDS ${MAKEDEPS_INCLUDE}/.stamp)

# Build a Makefile-based dep against MAKEDEPS_PREFIX, then publish its public
# headers (and static lib) back into the same prefix so a downstream Makefile
# dep (e.g. gbwtgraph -> gbwt) can find them.
function(vg_makefile_dep target_name dep_name lib_relpath public_header_subdir)
    set(lib_src  "${VG_DEPS_DIR}/${dep_name}/${lib_relpath}")
    set(lib_dst  "${MAKEDEPS_LIB}/lib${dep_name}.a")
    set(hdr_src  "${VG_DEPS_DIR}/${dep_name}/include/${public_header_subdir}")
    set(hdr_dst  "${MAKEDEPS_INCLUDE}/${public_header_subdir}")

    ExternalProject_Add(${dep_name}_ep
        SOURCE_DIR        ${VG_DEPS_DIR}/${dep_name}
        BUILD_IN_SOURCE   1
        CONFIGURE_COMMAND ${CMAKE_COMMAND} -E make_directory
                              ${VG_DEPS_DIR}/${dep_name}/obj
            COMMAND       ${CMAKE_COMMAND} -E make_directory
                              ${VG_DEPS_DIR}/${dep_name}/lib
        BUILD_COMMAND     $(MAKE) ${lib_relpath}
            SDSL_DIR=${sdsl-lite_SOURCE_DIR}
            INC_DIR=${MAKEDEPS_INCLUDE}
            LIB_DIR=${MAKEDEPS_LIB}
        # Publish headers + lib into the staging prefix.
        INSTALL_COMMAND   ${CMAKE_COMMAND} -E copy_directory ${hdr_src} ${hdr_dst}
            COMMAND       ${CMAKE_COMMAND} -E copy ${lib_src} ${lib_dst}
        BUILD_BYPRODUCTS  ${lib_src} ${lib_dst}
        DEPENDS           makedeps_stage)

    add_library(${target_name} STATIC IMPORTED)
    set_target_properties(${target_name} PROPERTIES
        IMPORTED_LOCATION             ${lib_dst}
        INTERFACE_INCLUDE_DIRECTORIES "${MAKEDEPS_INCLUDE}")
    add_dependencies(${target_name} ${dep_name}_ep)
endfunction()

vg_makefile_dep(vg_gbwt      gbwt      lib/libgbwt.a      gbwt)
vg_makefile_dep(vg_gcsa2     gcsa2     lib/libgcsa2.a     gcsa)
# gbwtgraph #includes <gbwt/...> — must run after gbwt's headers are staged.
vg_makefile_dep(vg_gbwtgraph gbwtgraph lib/libgbwtgraph.a gbwtgraph)
add_dependencies(gbwtgraph_ep gbwt_ep)

# sonLib — C utility library. Headers stage into sonLib/ subdir for <sonLib/X.h> includes.
ExternalProject_Add(sonLib_ep
    SOURCE_DIR        ${VG_DEPS_DIR}/sonLib
    BUILD_IN_SOURCE   1
    CONFIGURE_COMMAND ${CMAKE_COMMAND} -E make_directory ${VG_DEPS_DIR}/sonLib/lib
              COMMAND ${CMAKE_COMMAND} -E make_directory ${VG_DEPS_DIR}/sonLib/C/obj
    BUILD_COMMAND     $(MAKE) -C C ../lib/sonLib.a
    INSTALL_COMMAND   ${CMAKE_COMMAND} -E make_directory ${MAKEDEPS_INCLUDE}/sonLib
              COMMAND ${CMAKE_COMMAND} -E copy_directory
                      ${VG_DEPS_DIR}/sonLib/C/inc ${MAKEDEPS_INCLUDE}/sonLib
              COMMAND ${CMAKE_COMMAND} -E copy_directory
                      ${VG_DEPS_DIR}/sonLib/lib/ ${MAKEDEPS_INCLUDE}/sonLib
              COMMAND ${CMAKE_COMMAND} -E copy
                      ${VG_DEPS_DIR}/sonLib/lib/sonLib.a ${MAKEDEPS_LIB}/libsonLib.a
    BUILD_BYPRODUCTS  ${VG_DEPS_DIR}/sonLib/lib/sonLib.a
                      ${MAKEDEPS_LIB}/libsonLib.a)

add_library(vg_sonlib STATIC IMPORTED)
set_target_properties(vg_sonlib PROPERTIES
    IMPORTED_LOCATION             ${MAKEDEPS_LIB}/libsonLib.a
    INTERFACE_INCLUDE_DIRECTORIES ${MAKEDEPS_INCLUDE})
add_dependencies(vg_sonlib sonLib_ep)

# jemalloc — only when selected. Built via its autotools.
if(VG_USE_JEMALLOC)
    set(JEMALLOC_INSTALL_DIR "${VG_STAGE_DIR}/jemalloc")
    ExternalProject_Add(jemalloc_ep
        SOURCE_DIR        ${VG_DEPS_DIR}/jemalloc/malloc
        BUILD_IN_SOURCE   1
        CONFIGURE_COMMAND ./autogen.sh --prefix=${JEMALLOC_INSTALL_DIR} --with-jemalloc-prefix=
        BUILD_COMMAND     $(MAKE)
        INSTALL_COMMAND   $(MAKE) install
        BUILD_BYPRODUCTS  ${JEMALLOC_INSTALL_DIR}/lib/libjemalloc.a)
    file(MAKE_DIRECTORY ${JEMALLOC_INSTALL_DIR}/include)
    add_library(vg_allocator STATIC IMPORTED GLOBAL)
    set_target_properties(vg_allocator PROPERTIES
        IMPORTED_LOCATION             ${JEMALLOC_INSTALL_DIR}/lib/libjemalloc.a
        INTERFACE_INCLUDE_DIRECTORIES ${JEMALLOC_INSTALL_DIR}/include)
    add_dependencies(vg_allocator jemalloc_ep)
    set(VG_ALLOCATOR_SOURCE src/config/allocator_config_jemalloc.cpp)
elseif(VG_USE_MIMALLOC)
    # mimalloc has a real CMakeLists.txt — pull it in.
    add_subdirectory(${VG_DEPS_DIR}/mimalloc ${CMAKE_BINARY_DIR}/deps/mimalloc EXCLUDE_FROM_ALL)
    add_library(vg_allocator INTERFACE)
    target_link_libraries(vg_allocator INTERFACE mimalloc-static)
    set(VG_ALLOCATOR_SOURCE src/config/allocator_config_mimalloc.cpp)
else()
    add_library(vg_allocator INTERFACE)
    set(VG_ALLOCATOR_SOURCE src/config/allocator_config_system.cpp)
endif()

# ---------------------------------------------------------------------------
# 3. Generated headers: stamp git version + build environment macros.
#    The hand-written src/version.cpp #includes these two headers and reads
#    VG_GIT_VERSION at compile time. Replaces the existing Makefile rule that
#    sed-generated the same files.
# ---------------------------------------------------------------------------
set(VG_GENERATED_DIR "${CMAKE_BINARY_DIR}/generated")
file(MAKE_DIRECTORY ${VG_GENERATED_DIR})

execute_process(
    COMMAND git describe --always --tags
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE VG_GIT_DESCRIBE
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET)
if(NOT VG_GIT_DESCRIBE)
    set(VG_GIT_DESCRIBE "unknown")
endif()

file(WRITE ${VG_GENERATED_DIR}/vg_git_version.hpp
    "#pragma once\n#define VG_GIT_VERSION \"${VG_GIT_DESCRIBE}\"\n")
file(WRITE ${VG_GENERATED_DIR}/vg_environment_version.hpp
    "#pragma once\n"
    "#define VG_COMPILER_VERSION \"${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}\"\n"
    "#define VG_OS \"${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_VERSION}\"\n"
    "#define VG_BUILD_USER \"$ENV{USER}\"\n"
    "#define VG_BUILD_HOST \"$ENV{HOSTNAME}\"\n")

# ---------------------------------------------------------------------------
# 4. vg_core: aggregate library of everything in src/ except main + subcommand
# ---------------------------------------------------------------------------
file(GLOB VG_CORE_TOPLEVEL CONFIGURE_DEPENDS src/*.cpp)
list(FILTER VG_CORE_TOPLEVEL EXCLUDE REGEX "/main\\.cpp$")

file(GLOB VG_CORE_ALGORITHMS CONFIGURE_DEPENDS src/algorithms/*.cpp)
file(GLOB VG_CORE_IO         CONFIGURE_DEPENDS src/io/*.cpp)

add_library(vg_core STATIC
    ${VG_CORE_TOPLEVEL}
    ${VG_CORE_ALGORITHMS}
    ${VG_CORE_IO}
    ${VG_ALLOCATOR_SOURCE})

target_include_directories(vg_core PUBLIC
    ${CMAKE_SOURCE_DIR}/src
    ${VG_GENERATED_DIR})

target_link_libraries(vg_core PUBLIC
    # vgteam (FetchContent)
    handlegraph_shared vgio_static libbdsg libgssw
    vg_xg progress_bar vcflib structures_static
    dozeu
    # vendored static libs
    vg_ssw vg_sha1 vg_fermi_lite vg_libdeflate vg_tabixpp vg_libvcfh
    vg_headers_only backward-cpp raptor2
    # external (ExternalProject)
    vg_htslib vg_gbwt vg_gcsa2 vg_gbwtgraph vg_sonlib
    # CMake-native locals
    sdsl
    # system
    OpenMP::OpenMP_CXX Threads::Threads
    ZLIB::ZLIB BZip2::BZip2 LibLZMA::LibLZMA
    vg_allocator)

if(TARGET vg_sublinear)
    target_link_libraries(vg_core PUBLIC vg_sublinear)
endif()
if(VG_USE_BACKWARD AND TARGET backward)
    target_link_libraries(vg_core PUBLIC backward)
endif()

# ---------------------------------------------------------------------------
# 5. Subcommand support + per-subcommand OBJECT libraries
# ---------------------------------------------------------------------------
add_library(vg_subcommand_support STATIC
    src/subcommand/subcommand.cpp
    src/subcommand/options.cpp)
target_link_libraries(vg_subcommand_support PUBLIC vg_core)

# Each *_main.cpp becomes its own OBJECT library. OBJECT libraries are linked
# TU-by-TU (not as archive members), so the static-constructor side effect
# that registers the subcommand survives the final link without any
# --whole-archive plumbing.
file(GLOB VG_SUBCOMMAND_MAINS CONFIGURE_DEPENDS src/subcommand/*_main.cpp)
set(VG_SUBCOMMAND_TARGETS "")
foreach(src IN LISTS VG_SUBCOMMAND_MAINS)
    get_filename_component(stem ${src} NAME_WE)            # e.g. "giraffe_main"
    string(REGEX REPLACE "_main$" "" name ${stem})         # -> "giraffe"
    set(tgt "vg_sub_${name}")
    add_library(${tgt} OBJECT ${src})
    target_link_libraries(${tgt} PUBLIC vg_subcommand_support)
    list(APPEND VG_SUBCOMMAND_TARGETS ${tgt})
endforeach()

# ---------------------------------------------------------------------------
# 6. Unit-test object library (the `vg test` subcommand needs these linked in)
# ---------------------------------------------------------------------------
if(VG_BUILD_UNITTESTS)
    file(GLOB VG_UNITTEST_SOURCES         CONFIGURE_DEPENDS src/unittest/*.cpp)
    file(GLOB VG_UNITTEST_SUPPORT_SOURCES CONFIGURE_DEPENDS src/unittest/support/*.cpp)
    add_library(vg_unittest OBJECT
        ${VG_UNITTEST_SOURCES}
        ${VG_UNITTEST_SUPPORT_SOURCES})
    target_link_libraries(vg_unittest PUBLIC vg_subcommand_support)
endif()

# ---------------------------------------------------------------------------
# 7. Final vg executable
# ---------------------------------------------------------------------------
add_executable(vg src/main.cpp)
target_link_libraries(vg PRIVATE
    ${VG_SUBCOMMAND_TARGETS}
    vg_subcommand_support
    vg_core)
if(VG_BUILD_UNITTESTS)
    target_link_libraries(vg PRIVATE vg_unittest)
endif()

# ---------------------------------------------------------------------------
# 8. Install + test wiring
# ---------------------------------------------------------------------------
install(TARGETS vg RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

enable_testing()
add_test(NAME unit COMMAND vg test)

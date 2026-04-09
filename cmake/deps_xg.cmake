# Phase 6: xg — compiled directly (bypasses xg's own CMake build system)
#
# The vg Makefile compiles xg with a single hand-written command:
#   $(CXX) $(INCLUDE_FLAGS) $(CXXFLAGS) -fPIC -DNO_GFAKLUGE -c -o xg.o xg.cpp
#   ar rs libxg.a xg.o
#
# We replicate this with ExternalProject_Add so we can express the dependency
# on all the header-producing deps (libhandlegraph, sdsl, mmmultimap, ips4o, mio).
# libbdsg is placed in Phase 5 (deps_external_project.cmake) since it has its
# own CMakeLists.txt. Only xg remains here.

include(ExternalProject)

set(DEPS_DIR ${CMAKE_SOURCE_DIR}/deps)
set(XG_DIR   ${DEPS_DIR}/xg)

ExternalProject_Add(xg_ep
    SOURCE_DIR  ${XG_DIR}
    BUILD_IN_SOURCE ON
    CONFIGURE_COMMAND
        # Copy xg's headers into the staging include dir
        ${CMAKE_COMMAND} -E copy_directory ${XG_DIR}/src ${VG_INC_DIR}
    BUILD_COMMAND
        ${CMAKE_CXX_COMPILER}
            -std=c++17
            -O3 -ggdb -g
            -fPIC
            -DNO_GFAKLUGE
            -I${VG_INC_DIR}
            -I${CMAKE_SOURCE_DIR}/src
            -c -o ${XG_DIR}/xg.o ${XG_DIR}/src/xg.cpp
        COMMAND ar rs ${VG_LIB_DIR}/libxg.a ${XG_DIR}/xg.o
    INSTALL_COMMAND ""
    # Must run after all header-producing deps are done
    DEPENDS
        libhandlegraph_ep
        sdsl_ep
        sparsehash_ep
        # Phase 3 header-only deps don't have ExternalProject targets;
        # file(COPY) runs at configure time so they're always available.
    BUILD_BYPRODUCTS
        ${VG_LIB_DIR}/libxg.a
)

add_library(dep_xg STATIC IMPORTED GLOBAL)
set_target_properties(dep_xg PROPERTIES
    IMPORTED_LOCATION             ${VG_LIB_DIR}/libxg.a
    INTERFACE_INCLUDE_DIRECTORIES ${VG_INC_DIR}
)
add_dependencies(dep_xg xg_ep)

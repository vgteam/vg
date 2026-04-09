# Phase 11: Utility targets and install rules
#
# Mirrors Makefile targets: lint, test, docs, man, static, clean targets.
# get-deps and set-path are not CMake concerns (system-level or env setup).

# ── lint ──────────────────────────────────────────────────────────────────
find_program(Python3_EXECUTABLE python3)

add_custom_target(lint
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/scripts/check_options.py
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    COMMENT "Running vg option consistency check (check_options.py)"
    VERBATIM
)

# ── docs ──────────────────────────────────────────────────────────────────
find_program(DOXYGEN_EXECUTABLE doxygen)

if(DOXYGEN_EXECUTABLE)
    add_custom_target(docs
        COMMAND ${DOXYGEN_EXECUTABLE}
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        COMMENT "Generating Doxygen documentation"
        VERBATIM
    )
else()
    add_custom_target(docs
        COMMAND ${CMAKE_COMMAND} -E echo "doxygen not found; skipping docs"
        COMMENT "doxygen not available"
    )
endif()

# ── man ───────────────────────────────────────────────────────────────────
# man page build: vgmanmd.py → doc/man/vg-manpage.md → pandoc → doc/man/vg.1
find_program(PANDOC_EXECUTABLE pandoc)

add_custom_command(
    OUTPUT ${CMAKE_SOURCE_DIR}/doc/man/vg-manpage.md
    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_SOURCE_DIR}/doc/man
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/doc/vgmanmd.py
            > ${CMAKE_SOURCE_DIR}/doc/man/vg-manpage.md
    DEPENDS vg ${CMAKE_SOURCE_DIR}/doc/vgmanmd.desc.md ${CMAKE_SOURCE_DIR}/doc/vgmanmd.py
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    COMMENT "Generating vg man page markdown"
    VERBATIM
)

add_custom_command(
    OUTPUT ${CMAKE_SOURCE_DIR}/doc/wiki/vg-manpage.md
    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_SOURCE_DIR}/doc/wiki
    # sed 1d: remove first line (wiki version strips the man-page-specific header)
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_SOURCE_DIR}/cmake/sed1d.cmake
    DEPENDS ${CMAKE_SOURCE_DIR}/doc/man/vg-manpage.md
    COMMENT "Generating wiki manpage markdown"
)
file(WRITE ${CMAKE_SOURCE_DIR}/cmake/sed1d.cmake
"file(READ \"${CMAKE_SOURCE_DIR}/doc/man/vg-manpage.md\" _content)
string(REGEX REPLACE \"^[^\\n]*\\n\" \"\" _stripped \"\${_content}\")
file(WRITE \"${CMAKE_SOURCE_DIR}/doc/wiki/vg-manpage.md\" \"\${_stripped}\")
")

if(PANDOC_EXECUTABLE)
    add_custom_command(
        OUTPUT ${CMAKE_SOURCE_DIR}/doc/man/vg.1
        COMMAND ${PANDOC_EXECUTABLE} --standalone --to man
                ${CMAKE_SOURCE_DIR}/doc/man/vg-manpage.md
                -o ${CMAKE_SOURCE_DIR}/doc/man/vg.1
        DEPENDS ${CMAKE_SOURCE_DIR}/doc/man/vg-manpage.md
        COMMENT "Converting man page markdown to roff"
        VERBATIM
    )
    add_custom_target(man
        DEPENDS
            ${CMAKE_SOURCE_DIR}/doc/man/vg.1
            ${CMAKE_SOURCE_DIR}/doc/wiki/vg-manpage.md
    )
else()
    add_custom_target(man
        COMMAND ${CMAKE_COMMAND} -E echo "pandoc not found; cannot build man page"
    )
endif()

# ── test (integration) ────────────────────────────────────────────────────
# Per-suite unit test binaries + CTest are set up in Phase 9.
# This 'integration_test' target mirrors the Makefile's 'make test' which runs
# the Perl prove suite and doc tests.
find_program(PROVE_EXECUTABLE prove)

add_custom_target(integration_test
    COMMAND ${CMAKE_COMMAND} -E echo "Running prove test suite..."
    COMMAND ${PROVE_EXECUTABLE} -v t
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/test
    DEPENDS vg
    COMMENT "Running vg integration tests (prove -v t)"
)

# doc tests (strip compiler env vars to avoid picking up vg's own libraries)
add_custom_target(doc_test
    COMMAND ${CMAKE_COMMAND} -E env
        CFLAGS= CXXFLAGS= CPPFLAGS= LDFLAGS= INCLUDE_FLAGS=
        LIBRARY_PATH= LD_LIBRARY_PATH= DYLD_LIBRARY_PATH=
        DYLD_FALLBACK_LIBRARY_PATH= LD_INCLUDE_PATH=
        CC= CXX= CXX_STANDARD=
        "PATH=${VG_BIN_DIR}:$ENV{PATH}"
        ${CMAKE_SOURCE_DIR}/doc/test-docs.sh
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    DEPENDS vg
    COMMENT "Running doc tests"
)

# ── test/build_graph helper binary ────────────────────────────────────────
# A small test binary built from test/build_graph.cpp
if(EXISTS ${CMAKE_SOURCE_DIR}/test/build_graph.cpp)
    add_executable(build_graph ${CMAKE_SOURCE_DIR}/test/build_graph.cpp)
    set_target_properties(build_graph PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/test
    )
    target_include_directories(build_graph PRIVATE ${VG_INC_DIR} ${VG_SRC_DIR})
    target_link_libraries(build_graph PRIVATE libvg vg_system_libs)
    # Makefile strips debug symbols due to macOS antimalware (HX Antimalware)
    if(APPLE)
        add_custom_command(TARGET build_graph POST_BUILD
            COMMAND strip ${CMAKE_SOURCE_DIR}/test/build_graph
            COMMENT "Stripping build_graph (macOS antimalware workaround)"
        )
    endif()
endif()

# ── Static build variant ──────────────────────────────────────────────────
# The Makefile's 'make static' re-links the vg binary with:
#   -static -static-libstdc++ -static-libgcc -Wl,--allow-multiple-definition
# We implement this as a separate 'vg_static' executable target
# with EXCLUDE_FROM_ALL so it doesn't build by default.
if(NOT APPLE)
    # Static builds are only supported on Linux (Makefile comment)
    add_executable(vg_static EXCLUDE_FROM_ALL
        ${VG_SRC_DIR}/main.cpp
        ${SUBCOMMAND_SRCS}
        ${UNITTEST_SRCS}
        ${UNITTEST_SUPPORT_SRCS}
        $<TARGET_OBJECTS:allocator_config>
    )
    set_target_properties(vg_static PROPERTIES
        OUTPUT_NAME vg
        RUNTIME_OUTPUT_DIRECTORY ${VG_BIN_DIR}
    )
    target_include_directories(vg_static PRIVATE
        ${VG_INC_DIR} ${VG_INC_DIR}/dynamic
        ${CMAKE_SOURCE_DIR} ${VG_SRC_DIR}
        ${VG_SRC_DIR}/unittest ${VG_SRC_DIR}/unittest/support
        ${VG_SRC_DIR}/subcommand
    )
    target_link_options(vg_static PRIVATE
        -static
        -static-libstdc++
        -static-libgcc
        -Wl,--allow-multiple-definition
        -L${VG_LIB_DIR}
    )
    # Same libraries as the dynamic vg, but all static
    target_link_libraries(vg_static PRIVATE
        libvg
        dep_vcflib dep_wfa2 dep_tabixpp dep_gssw dep_ssw dep_sublinearLS
        pthread ncurses
        dep_gcsa2 dep_gbwtgraph dep_gbwt dep_kff
        dep_divsufsort dep_divsufsort64 dep_libvcfh dep_raptor
        dep_pinchesandcacti dep_3edgeconnected dep_sonlib dep_fml dep_structures
        dep_libbdsg dep_xg dep_sdsl dep_libhandlegraph crypto
        dep_dwfl dep_dw dep_dwelf dep_ebl dep_elf atomic
        ${VG_OMP_LINK_LIBS} Boost::program_options PkgConfig::CAIRO PkgConfig::ZSTD
        dep_libvgio dep_htslib dep_libdeflate z bz2 lzma
        protobuf::libprotobuf PkgConfig::JANSSON
        pthread m
        $<$<BOOL:${VG_USE_JEMALLOC}>:${VG_LIB_DIR}/libjemalloc.a>
    )

    add_custom_target(static DEPENDS vg_static
        COMMENT "Building static vg binary"
    )

    add_custom_target(static_docker
        COMMAND strip -d ${VG_BIN_DIR}/vg
        COMMAND ${CMAKE_COMMAND} -E env DOCKER_BUILDKIT=1
            docker build ${CMAKE_SOURCE_DIR} -f ${CMAKE_SOURCE_DIR}/Dockerfile.static -t vg
        DEPENDS static
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        COMMENT "Building static Docker image"
    )
endif()

# ── shuf symlink (test helper) ────────────────────────────────────────────
# Makefile: ln -s `which gshuf` bin/shuf  (Mac) or  ln -s `which shuf` bin/shuf  (Linux)
if(APPLE)
    find_program(_shuf_src gshuf)
else()
    find_program(_shuf_src shuf)
endif()
if(_shuf_src)
    add_custom_target(shuf_symlink
        COMMAND ${CMAKE_COMMAND} -E create_symlink ${_shuf_src} ${VG_BIN_DIR}/shuf
        COMMENT "Creating shuf symlink in bin/"
    )
endif()

# ── clean targets (for reference) ─────────────────────────────────────────
# CMake's built-in 'cmake --build . --target clean' removes all build artifacts.
# The Makefile's 'clean-vg' (vg objects only) and 'clean-vcflib' (vcflib only)
# have no direct CMake equivalent — use the ExternalProject clean targets or
# remove specific directories from the build tree manually.

# ── install rules ─────────────────────────────────────────────────────────
include(GNUInstallDirs)

install(TARGETS vg
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
install(TARGETS libvg libvg_shared
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
install(DIRECTORY ${VG_INC_DIR}/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# Phase 9: Per-suite unit test binaries
#
# Makefile pattern:
#   $(UNITTEST_EXE): $(UNITTEST_BIN_DIR)/%: ...
#       $(CXX) ... -o $@ $<
#           $(UNITTEST_SUPPORT_OBJ) $(CONFIG_OBJ)
#           $(LIB_DIR)/libvg.$(SHARED_SUFFIX)
#           $(LD_LIB_FLAGS) $(LD_STATIC_LIB_FLAGS) $(LD_STATIC_LIB_DEPS)
#
# Each test suite .cpp gets its own small binary.  They all link against the
# SHARED libvg (not the static one) so individual test builds are fast.
# driver.cpp (which has main()) IS included here — unlike the main vg binary.

set(VG_SRC_DIR ${CMAKE_SOURCE_DIR}/src)

# All support sources including driver.cpp (provides main() for test binaries)
file(GLOB UNITTEST_SUPPORT_ALL_SRCS
    ${VG_SRC_DIR}/unittest/support/*.cpp
)

# Determine the set of test suites (one binary per .cpp in src/unittest/)
file(GLOB UNITTEST_CPPS ${VG_SRC_DIR}/unittest/*.cpp)

enable_testing()

foreach(test_src IN LISTS UNITTEST_CPPS)
    get_filename_component(test_name ${test_src} NAME_WE)
    set(exe_name vg_test_${test_name})

    add_executable(${exe_name}
        ${test_src}
        ${UNITTEST_SUPPORT_ALL_SRCS}
        $<TARGET_OBJECTS:allocator_config>
    )

    set_target_properties(${exe_name} PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${VG_BIN_DIR}/unittest
    )

    target_include_directories(${exe_name} PRIVATE
        ${VG_INC_DIR}
        ${VG_INC_DIR}/dynamic
        ${CMAKE_SOURCE_DIR}
        ${VG_SRC_DIR}
        ${VG_SRC_DIR}/unittest
        ${VG_SRC_DIR}/unittest/support
        ${VG_SRC_DIR}/subcommand
    )

    # Link against the shared libvg for fast per-test builds (mirrors Makefile)
    target_link_libraries(${exe_name} PRIVATE
        libvg_shared

        # LD_LIB_FLAGS (no Bstatic wrapper for test binaries)
        dep_vcflib dep_wfa2 dep_tabixpp dep_gssw dep_ssw dep_sublinearLS
        pthread ncurses
        dep_gcsa2 dep_gbwtgraph dep_gbwt dep_kff
        dep_divsufsort dep_divsufsort64
        dep_libvcfh dep_raptor dep_pinchesandcacti dep_3edgeconnected
        dep_sonlib dep_fml dep_structures dep_libbdsg dep_xg dep_sdsl
        dep_libhandlegraph
        crypto

        # LD_STATIC_LIB_FLAGS (no Bstatic wrapper — test binaries are dynamic)
        dep_libvgio dep_htslib dep_libdeflate z bz2 lzma
        protobuf::libprotobuf PkgConfig::JANSSON

        # LD_STATIC_LIB_DEPS
        pthread m

        # LD_EXE_LIB_FLAGS
        $<$<BOOL:${VG_USE_JEMALLOC}>:${VG_LIB_DIR}/libjemalloc.a>

        # System libs
        PkgConfig::CAIRO PkgConfig::ZSTD
        Boost::program_options
        ${VG_OMP_LINK_LIBS}
    )

    if(NOT APPLE)
        target_link_libraries(${exe_name} PRIVATE
            dep_dwfl dep_dw dep_dwelf dep_ebl dep_elf
            atomic
        )
    endif()

    target_link_options(${exe_name} PRIVATE
        -L${VG_LIB_DIR}
        LINKER:-rpath,${VG_LIB_DIR}
    )

    # Register with CTest (runs the binary with no args = run all tests in suite)
    add_test(NAME ${test_name} COMMAND ${exe_name})
    set_tests_properties(${test_name} PROPERTIES
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/test
    )
endforeach()

# Convenience target: build all test suite binaries
add_custom_target(vg_test_all
    DEPENDS ${UNITTEST_CPPS}
    COMMENT "Building all per-suite vg unit test binaries"
)
foreach(test_src IN LISTS UNITTEST_CPPS)
    get_filename_component(test_name ${test_src} NAME_WE)
    add_dependencies(vg_test_all vg_test_${test_name})
endforeach()

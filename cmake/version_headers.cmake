# Phase 10: Version header generation
#
# Mirrors the Makefile's global $(shell ...) calls that generate:
#   src/vg_git_version.hpp
#   src/vg_environment_version.hpp
#
# These are generated at CMake configure time (equivalent to the Makefile
# evaluating them at parse time with $(shell ...)).
#
# The "only update if changed" behavior (Makefile's diff trick) is replicated
# by comparing the new content to the existing file and only writing if different.

set(VG_SRC_DIR ${CMAKE_SOURCE_DIR}/src)

# ── Git version ──────────────────────────────────────────────────────────
execute_process(
    COMMAND git describe --always --tags
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE _vg_git_version
    ERROR_VARIABLE  _vg_git_error
    OUTPUT_STRIP_TRAILING_WHITESPACE
    RESULT_VARIABLE _vg_git_result
)

if(NOT _vg_git_result EQUAL 0 OR _vg_git_version STREQUAL "")
    set(_vg_git_version "git-error")
    message(WARNING "Could not determine git version: ${_vg_git_error}")
endif()

set(_git_content "#define VG_GIT_VERSION \"${_vg_git_version}\"\n")

# Write only if content changed (mirrors Makefile's diff trick)
set(_git_header ${VG_SRC_DIR}/vg_git_version.hpp)
if(EXISTS ${_git_header})
    file(READ ${_git_header} _existing_git)
else()
    set(_existing_git "")
endif()

if(NOT _existing_git STREQUAL _git_content)
    file(WRITE ${_git_header} ${_git_content})
    message(STATUS "Updated vg_git_version.hpp: ${_vg_git_version}")
else()
    message(STATUS "vg_git_version.hpp unchanged: ${_vg_git_version}")
endif()

# ── Environment version ───────────────────────────────────────────────────
execute_process(
    COMMAND ${CMAKE_CXX_COMPILER} --version
    OUTPUT_VARIABLE _cxx_version_full
    ERROR_VARIABLE  _cxx_version_err
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
# Take just the first line (mirrors Makefile's | head -n 1)
string(REGEX REPLACE "\n.*$" "" _cxx_version "${_cxx_version_full}")

execute_process(
    COMMAND uname
    OUTPUT_VARIABLE _vg_os
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

execute_process(
    COMMAND whoami
    OUTPUT_VARIABLE _vg_build_user
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

execute_process(
    COMMAND hostname
    OUTPUT_VARIABLE _vg_build_host
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

set(_env_content
"#define VG_COMPILER_VERSION \"${_cxx_version}\"
#define VG_OS \"${_vg_os}\"
#define VG_BUILD_USER \"${_vg_build_user}\"
#define VG_BUILD_HOST \"${_vg_build_host}\"
")

set(_env_header ${VG_SRC_DIR}/vg_environment_version.hpp)
if(EXISTS ${_env_header})
    file(READ ${_env_header} _existing_env)
else()
    set(_existing_env "")
endif()

if(NOT _existing_env STREQUAL _env_content)
    file(WRITE ${_env_header} ${_env_content})
    message(STATUS "Updated vg_environment_version.hpp")
else()
    message(STATUS "vg_environment_version.hpp unchanged")
endif()

# ── Tell CMake that version.cpp depends on these generated headers ─────────
# (CMake tracks .cpp → .hpp dependencies automatically if the .hpp is
#  included via target_include_directories, but we make it explicit here
#  for documentation and correctness.)
set_source_files_properties(
    ${VG_SRC_DIR}/version.cpp
    PROPERTIES OBJECT_DEPENDS
        "${VG_SRC_DIR}/vg_git_version.hpp;${VG_SRC_DIR}/vg_environment_version.hpp"
)

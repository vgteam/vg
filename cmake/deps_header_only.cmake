# Phase 3: Header-only and micro-compile dependencies
#
# For deps where the Makefile simply copies headers, we replicate that here at
# CMake configure time using file(COPY ...).  Each dep also gets an INTERFACE
# target so downstream targets can express the dependency cleanly.
#
# Deps covered:
#   lru_cache, sparsepp, mio, atomic_queue, backward-cpp, dozeu,
#   ips4o, BBHash, mmmultimap, DYNAMIC (+ hopscotch_map)
#   sha1, progress_bar  (compiled into libvg; see Phase 7)
#
# NOTE: sparsehash needs ./configure (generates sparseconfig.h).
#       It is handled in Phase 5 (autoconf deps) as dep_sparsehash.

set(DEPS_DIR ${CMAKE_SOURCE_DIR}/deps)

# ── lru_cache ──────────────────────────────────────────────────────────────
# Makefile: cp deps/lru_cache/*.h* include/
file(GLOB _lru_hdrs ${DEPS_DIR}/lru_cache/*.h ${DEPS_DIR}/lru_cache/*.h)
file(COPY ${_lru_hdrs} DESTINATION ${VG_INC_DIR})

add_library(dep_lru_cache INTERFACE)
target_include_directories(dep_lru_cache INTERFACE ${VG_INC_DIR})

# ── sparsepp ───────────────────────────────────────────────────────────────
# Makefile: cp -r deps/sparsepp/sparsepp include/
file(COPY ${DEPS_DIR}/sparsepp/sparsepp DESTINATION ${VG_INC_DIR})

add_library(dep_sparsepp INTERFACE)
target_include_directories(dep_sparsepp INTERFACE ${VG_INC_DIR})

# ── mio ────────────────────────────────────────────────────────────────────
# Makefile: cp -r deps/mio/include/mio include/mio/
file(COPY ${DEPS_DIR}/mio/include/mio DESTINATION ${VG_INC_DIR})

add_library(dep_mio INTERFACE)
target_include_directories(dep_mio INTERFACE ${VG_INC_DIR})

# ── atomic_queue ───────────────────────────────────────────────────────────
# Makefile: cp -r deps/atomic_queue/include/atomic_queue/* include/
# (flattens the directory into include/ so code does #include <atomic_queue.h>)
file(GLOB _aq_hdrs ${DEPS_DIR}/atomic_queue/include/atomic_queue/*)
file(COPY ${_aq_hdrs} DESTINATION ${VG_INC_DIR})

add_library(dep_atomic_queue INTERFACE)
target_include_directories(dep_atomic_queue INTERFACE ${VG_INC_DIR})

# ── backward-cpp ──────────────────────────────────────────────────────────
# Makefile: cp deps/backward-cpp/backward.hpp include/
file(COPY ${DEPS_DIR}/backward-cpp/backward.hpp DESTINATION ${VG_INC_DIR})

add_library(dep_backward_cpp INTERFACE)
target_include_directories(dep_backward_cpp INTERFACE ${VG_INC_DIR})

# ── dozeu ──────────────────────────────────────────────────────────────────
# Makefile: cp -r deps/dozeu/simde include/simde
#           mkdir include/dozeu && cp deps/dozeu/*.h include/dozeu/
file(COPY ${DEPS_DIR}/dozeu/simde DESTINATION ${VG_INC_DIR})
file(MAKE_DIRECTORY ${VG_INC_DIR}/dozeu)
file(GLOB_RECURSE _dozeu_hdrs ${DEPS_DIR}/dozeu/*.h)
file(COPY ${_dozeu_hdrs} DESTINATION ${VG_INC_DIR}/dozeu)

add_library(dep_dozeu INTERFACE)
target_include_directories(dep_dozeu INTERFACE ${VG_INC_DIR})

# ── ips4o ──────────────────────────────────────────────────────────────────
# Makefile: cp -r deps/ips4o/ips4o* include/
file(GLOB _ips4o_items ${DEPS_DIR}/ips4o/ips4o.hpp ${DEPS_DIR}/ips4o/ips4o)
foreach(_item IN LISTS _ips4o_items)
    if(IS_DIRECTORY ${_item})
        file(COPY ${_item} DESTINATION ${VG_INC_DIR})
    else()
        file(COPY ${_item} DESTINATION ${VG_INC_DIR})
    endif()
endforeach()

add_library(dep_ips4o INTERFACE)
target_include_directories(dep_ips4o INTERFACE ${VG_INC_DIR})

# ── BBHash ─────────────────────────────────────────────────────────────────
# Makefile: cp deps/BBHash/BooPHF.h include/
file(COPY ${DEPS_DIR}/BBHash/BooPHF.h DESTINATION ${VG_INC_DIR})

add_library(dep_bbhash INTERFACE)
target_include_directories(dep_bbhash INTERFACE ${VG_INC_DIR})

# ── mmmultimap ─────────────────────────────────────────────────────────────
# Makefile: cp deps/mmmultimap/src/mmmultimap.hpp deps/mmmultimap/src/mmmultiset.hpp include/
file(COPY
    ${DEPS_DIR}/mmmultimap/src/mmmultimap.hpp
    ${DEPS_DIR}/mmmultimap/src/mmmultiset.hpp
    DESTINATION ${VG_INC_DIR}
)

add_library(dep_mmmultimap INTERFACE)
target_include_directories(dep_mmmultimap INTERFACE
    ${VG_INC_DIR}
    # mmmultimap pulls in mio and atomic_queue transitively
)
target_link_libraries(dep_mmmultimap INTERFACE dep_mio dep_atomic_queue)

# ── DYNAMIC (+ hopscotch_map) ──────────────────────────────────────────────
# Makefile runs cmake in deps/DYNAMIC to build, then copies:
#   deps/DYNAMIC/deps/hopscotch_map/include/* → include/
#   deps/DYNAMIC/include/dynamic/*            → include/dynamic/
#
# DYNAMIC itself is header-only. We point include paths directly at the
# submodule to avoid running cmake just for headers.
#
# Code includes:
#   #include "dynamic/dynamic.hpp"   → deps/DYNAMIC/include/dynamic/dynamic.hpp
#   #include <tsl/hopscotch_map.h>   → deps/DYNAMIC/deps/hopscotch_map/include/tsl/...

add_library(dep_dynamic INTERFACE)
target_include_directories(dep_dynamic INTERFACE
    ${DEPS_DIR}/DYNAMIC/include
    ${DEPS_DIR}/DYNAMIC/deps/hopscotch_map/include
)

# ── sha1 (compiled, linked into libvg) ────────────────────────────────────
# Makefile: obj/sha1.o from deps/sha1/sha1.cpp; header copied to include/
# Defined here as an OBJECT library; added to libvg in Phase 7.
file(COPY
    ${DEPS_DIR}/sha1/sha1.hpp
    DESTINATION ${VG_INC_DIR}
)

add_library(dep_sha1 OBJECT ${DEPS_DIR}/sha1/sha1.cpp)
target_include_directories(dep_sha1 PRIVATE ${VG_INC_DIR})

# ── progress_bar (compiled, linked into libvg) ────────────────────────────
# Makefile: obj/progress_bar.o from deps/progress_bar/progress_bar.cpp
file(COPY ${DEPS_DIR}/progress_bar/progress_bar.hpp DESTINATION ${VG_INC_DIR})

add_library(dep_progress_bar OBJECT ${DEPS_DIR}/progress_bar/progress_bar.cpp)
target_include_directories(dep_progress_bar PRIVATE ${VG_INC_DIR})

# ── Convenience aggregate target ──────────────────────────────────────────
# Any target can link dep_headers_all to get all header-only include paths.
add_library(dep_headers_all INTERFACE)
target_link_libraries(dep_headers_all INTERFACE
    dep_lru_cache
    dep_sparsepp
    dep_mio
    dep_atomic_queue
    dep_backward_cpp
    dep_dozeu
    dep_ips4o
    dep_bbhash
    dep_mmmultimap
    dep_dynamic
)

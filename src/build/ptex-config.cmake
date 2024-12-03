@PACKAGE_INIT@
include("${CMAKE_CURRENT_LIST_DIR}/ptex-exports.cmake")

include(CMakeFindDependencyMacro)

set(CMAKE_THREAD_PREFER_PTHREAD ON)
set(THREADS_PREFER_PTHREAD_FLAG ON)

find_dependency(Threads REQUIRED)

find_package(ZLIB REQUIRED)

set_and_check(Ptex_DIR @PACKAGE_CMAKE_INSTALL_PREFIX@)
set_and_check(Ptex_INCLUDE_DIRS @PACKAGE_CMAKE_INSTALL_INCLUDEDIR@)

# set(Ptex_FOUND TRUE)
set(PTEX_FOUND TRUE)

if( TARGET Ptex::Ptex_dynamic)
    set(PTEX_LIBRARY Ptex::Ptex_dynamic)
elseif(TARGET Ptex::Ptex_static)
    set(PTEX_LIBRARY Ptex::Ptex_static)
endif()

get_target_property(PTEX_INCLUDE_DIR ${PTEX_LIBRARY} INTERFACE_INCLUDE_DIRECTORIES)


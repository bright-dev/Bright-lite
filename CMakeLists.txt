CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
# taken from http://geant4.cern.ch/support/source/geant4/CMakeLists.txt
IF(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
    message(STATUS "Reactor requires an out-of-source build.")
    message(STATUS "Please remove these files from ${CMAKE_BINARY_DIR} first:")
    message(STATUS "CMakeCache.txt")
    message(STATUS "CMakeFiles")
    message(STATUS "Once these files are removed, create a separate directory")
    message(STATUS "and run CMake from there")
    message(FATAL_ERROR "in-source build detected")
ENDIF()

# This project name is reactor.
PROJECT( reactor )

# quiets fortify_source warnings when not compiling with optimizations
# in linux distros where compilers were compiled with fortify_source enabled by
# default (e.g. Arch linux).
MESSAGE("${CMAKE_BUILD_TYPE}")
STRING(TOLOWER "${CMAKE_BUILD_TYPE}" BUILD_TYPE)
IF(NOT ${BUILD_TYPE} STREQUAL "release")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -U_FORTIFY_SOURCE -D_FORTIFY_SOURCE=0")
ENDIF()
# no overflow warnings because of silly coin-ness
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-overflow")

# Direct any out-of-source builds to this directory
SET( REACTOR_SOURCE_DIR ${CMAKE_SOURCE_DIR} )

# Direct any binary installation paths to this directory
SET( REACTOR_BINARY_DIR ${CMAKE_BINARY_DIR} )

# This makes all the libraries build as SHARED
SET(BUILD_SHARED_LIBS true)

# Setup build locations.
IF(NOT CMAKE_RUNTIME_OUTPUT_DIRECTORY)
  SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${REACTOR_BINARY_DIR}/bin)
endif()
IF(NOT CMAKE_LIBRARY_OUTPUT_DIRECTORY)
  SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${REACTOR_BINARY_DIR}/lib)
endif()
IF(NOT CMAKE_ARCHIVE_OUTPUT_DIRECTORY)
  SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${REACTOR_BINARY_DIR}/lib)
ENDIF()

SET(REACTOR_EXECUTABLE_DIR ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

# use, i.e. don't skip the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib:${CMAKE_INSTALL_PREFIX}/lib/cyclus")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# the RPATH to be used when installing, but only if it's not a system directory
LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
IF("${isSystemDir}" STREQUAL "-1")
    SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib:${CMAKE_INSTALL_PREFIX}/lib/cyclus")
    GET_FILENAME_COMPONENT(cxxCompilerRoot ${CMAKE_CXX_COMPILER} DIRECTORY)
    GET_FILENAME_COMPONENT(cxxCompilerRoot ${cxxCompilerRoot} DIRECTORY)
    IF (NOT "${CMAKE_INSTALL_RPATH}" STREQUAL "${cxxCompilerRoot}")
        SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_RPATH}:${cxxCompilerRoot}/lib:${cxxCompilerRoot}/lib/cyclus")
    ENDIF (NOT "${CMAKE_INSTALL_RPATH}" STREQUAL "${cxxCompilerRoot}")
ENDIF("${isSystemDir}" STREQUAL "-1")
MESSAGE("-- CMAKE_INSTALL_RPATH: ${CMAKE_INSTALL_RPATH}")

# Tell CMake where the modules are
SET( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_DIR}/share/cmake-2.8/Modules")
SET( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${REACTOR_SOURCE_DIR}/cmake )
MESSAGE("CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}")


##############################################################################################
################################### begin find libraries #####################################
##############################################################################################

if (NOT DEPS_ROOT_DIR)
    get_filename_component(compdir ${CMAKE_CXX_COMPILER} DIRECTORY)
    get_filename_component(DEPS_ROOT_DIR ${compdir} DIRECTORY)
endif (NOT DEPS_ROOT_DIR)
SET(DEPS_HINTS HINTS "${DEPS_ROOT_DIR}")
SET(DEPS_BIN_HINTS ${DEPS_HINTS} "${DEPS_ROOT_DIR}/bin")
SET(DEPS_LIB_HINTS ${DEPS_HINTS} "${DEPS_ROOT_DIR}/lib")
SET(DEPS_INCLUDE_HINTS HINTS "${DEPS_ROOT_DIR}/include")
MESSAGE("-- Dependency Root Dir (DEPS_ROOT_DIR): ${DEPS_ROOT_DIR}")
MESSAGE("-- Dependency Hints (DEPS_HINTS): ${DEPS_HINTS}")
MESSAGE("-- Dependency Binary Hints (DEPS_BIN_HINTS): ${DEPS_BIN_HINTS}")
MESSAGE("-- Dependency Library Hints (DEPS_LIB_HINTS): ${DEPS_LIB_HINTS}")
MESSAGE("-- Dependency Include Hints (DEPS_INCLUDE_HINTS): ${DEPS_INCLUDE_HINTS}")

# Include macros
INCLUDE(UseCyclus)

# Find cyclus
FIND_PACKAGE( Cyclus REQUIRED)
SET(REACTOR_INCLUDE_DIRS ${REACTOR_INCLUDE_DIRS} ${CYCLUS_CORE_INCLUDE_DIR})
SET(LIBS ${LIBS} ${CYCLUS_CORE_LIBRARIES})

# Debian installs useful LibXML2 files to /usr/include/libxml2/libxml
# libxml2 is required for relaxng schema validation
FIND_PACKAGE(LibXml2 REQUIRED ${DEPS_HINTS})
ADD_DEFINITIONS(${LIBXML2_DEFINITIONS})
SET(LIBS ${LIBS} ${LIBXML2_LIBRARIES})

# Find LibXML++ and dependencies
FIND_PACKAGE( LibXML++ REQUIRED )
SET(REACTOR_INCLUDE_DIRS ${REACTOR_INCLUDE_DIRS} ${LibXML++_INCLUDE_DIR}
                         ${Glibmm_INCLUDE_DIRS} ${LibXML++Config_INCLUDE_DIR})
SET(LIBS ${LIBS} ${LibXML++_LIBRARIES})
MESSAGE("--LIBS: ${LIBS}")

MESSAGE("--LD_LIBRARY_PATH: $ENV{LD_LIBRARY_PATH}")

# Include the boost header files, system, and filesystem libraries
SET(Boost_USE_STATIC_LIBS       OFF)
SET(Boost_USE_STATIC_RUNTIME    OFF)
if (DEPS_ROOT_DIR)
  SET(BOOST_ROOT "${DEPS_ROOT_DIR}")
  SET(BOOST_INCLUDEDIR "${DEPS_INCLUDE_HINTS}")
endif (DEPS_ROOT_DIR)
FIND_PACKAGE( Boost COMPONENTS program_options filesystem system serialization REQUIRED)
SET(REACTOR_INCLUDE_DIRS ${REACTOR_INCLUDE_DIRS} ${Boost_INCLUDE_DIR})
MESSAGE("--    Boost Root: ${Boost_ROOT}")
MESSAGE("--    Boost Include directory: ${Boost_INCLUDE_DIR}")
MESSAGE("--    Boost Library directories: ${Boost_LIBRARY_DIRS}")
SET(LIBS ${LIBS} ${Boost_PROGRAM_OPTIONS_LIBRARY})
MESSAGE("--    Boost Program Options location: ${Boost_PROGRAM_OPTIONS_LIBRARY}")
SET(LIBS ${LIBS} ${Boost_SYSTEM_LIBRARY})
MESSAGE("--    Boost System location: ${Boost_SYSTEM_LIBRARY}")
SET(LIBS ${LIBS} ${Boost_FILESYSTEM_LIBRARY})
MESSAGE("--    Boost Filesystem location: ${Boost_FILESYSTEM_LIBRARY}")
SET(LIBS ${LIBS} ${Boost_SERIALIZATION_LIBRARY})
MESSAGE("--    Boost Serialization location: ${Boost_SERIALIZATION_LIBRARY}")

# find lapack and link to it
FIND_PACKAGE( LAPACK REQUIRED )
set(LIBS ${LIBS} ${LAPACK_LIBRARIES})
MESSAGE("-- Found LAPACK Linker Flags: ${LAPACK_LINKER_FLAGS}")
MESSAGE("-- Found LAPACK Libraries: ${LAPACK_LIBRARIES}")
MESSAGE("-- Found BLAS Libraries: ${BLAS_LIBRARIES}")

# Find HDF5
# Find HDF5
FIND_PACKAGE(HDF5 REQUIRED COMPONENTS HL)
ADD_DEFINITIONS(${HDF5_DEFINITIONS})
SET(REACTOR_INCLUDE_DIRS ${REACTOR_INCLUDE_DIRS} ${HDF5_INCLUDE_DIR})
set(LIBS ${LIBS} ${HDF5_LIBRARIES} ${HDF5_C_HL_LIBRARIES})
MESSAGE("--    HDF5 Root: ${HDF5_ROOT}")
MESSAGE("--    HDF5 Include directory: ${HDF5_INCLUDE_DIR}")
MESSAGE("--    HDF5 Library directories: ${HDF5_LIBRARY_DIRS}")
MESSAGE("--    HDF5 Libraries: ${HDF5_C_LIBRARIES}")
MESSAGE("--    HDF5 High Level Libraries: ${HDF5_C_HL_LIBRARIES}")

# find coin and link to it
FIND_PACKAGE( COIN REQUIRED )
SET(REACTOR_INCLUDE_DIRS ${REACTOR_INCLUDE_DIRS} ${COIN_INCLUDE_DIRS})
set(LIBS ${LIBS} ${COIN_LIBRARIES})
MESSAGE("--    COIN Root: ${COIN_ROOT}")
MESSAGE("--    COIN Include directories: ${COIN_INCLUDE_DIRS}")
MESSAGE("--    COIN Libraries: ${COIN_LIBRARIES}")

# include all the directories we just found
INCLUDE_DIRECTORIES(${REACTOR_INCLUDE_DIRS})


# add the agents
ADD_SUBDIRECTORY(src)

# uninstall target
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/cmake_uninstall.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
    IMMEDIATE @ONLY)

add_custom_target(uninstall
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)

# Install Brightlite Reactor Libraries
# finds reactor libraries based on presence of 'manifest.txt' file
set(BRIGHTLITE_SHARE "${CYCLUS_ROOT_DIR}/share/brightlite")
file(GLOB_RECURSE BLMANIFESTS FOLLOW_SYMLINKS "${CMAKE_SOURCE_DIR}/*/manifest.txt")
foreach(manifest ${BLMANIFESTS})
    get_filename_component(RXDIR ${manifest} DIRECTORY)
    set(RXDIRS "${RXDIRS}" "${RXDIR}")
endforeach()
install(DIRECTORY ${RXDIRS} DESTINATION "${BRIGHTLITE_SHARE}")

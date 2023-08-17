# Copyright 2023 TikTok Pte. Ltd.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

cmake_minimum_required(VERSION 3.14)

#######################################################
# Project DUET includes the following components:     #
#   1. DUET C++ library                               #
#   2. DUET C++ test                                  #
#   3. DUET C++ bench                                 #
#   4. DUET C++ example                               #
#######################################################

# [OPTION] CMAKE_BUILD_TYPE (DEFAULT: "Release")
# Select from Release, Debug, MiniSizeRel, or RelWithDebInfo.
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY
        STRINGS "Release" "Debug" "MinSizeRel" "RelWithDebInfo")
endif()
message(STATUS "Build type (CMAKE_BUILD_TYPE): ${CMAKE_BUILD_TYPE}")

project(DUET VERSION 0.1.0 LANGUAGES CXX C)

########################
# Global configuration #
########################

# CMake modules
include(CMakeDependentOption)
include(CMakePushCheckState)
include(CheckIncludeFiles)
include(CheckCXXSourceCompiles)
include(CheckCXXSourceRuns)

# Custom modules
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)
include(DuetCustomMacros)

# In Debug mode, define DUET_DEBUG.
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(DUET_DEBUG ON)
else()
    set(DUET_DEBUG OFF)
endif()
message(STATUS "DUET debug mode: ${DUET_DEBUG}")

# In Debug mode, enable extra compiler flags.
include(EnableDebugFlags)

# Always build position-independent-code
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# [OPTION] DUET_USE_CXX17 (default: ON)
# Use C++17, use C++14 otherwise.
set(DUET_USE_CXX17_OPTION_STR "Use C++17")
option(DUET_USE_CXX17 ${DUET_USE_CXX17_OPTION_STR} OFF)
message(STATUS "DUET_USE_CXX17: ${DUET_USE_CXX17}")
# Enable features from C++17 if available, disable features if set to OFF.
include(EnableCXX17)

# Add default files and directories.
include(GNUInstallDirs)

# Runtime path
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# Source Tree
set(DUET_INCLUDES_DIR ${CMAKE_CURRENT_LIST_DIR}/src)
set(DUET_CONFIG_IN_FILENAME ${CMAKE_CURRENT_LIST_DIR}/cmake/PETAce-DuetConfig.cmake.in)
set(DUET_CONFIG_H_IN_FILENAME ${DUET_INCLUDES_DIR}/duet/util/config.h.in)

# Build tree
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
set(DUET_CONFIG_FILENAME ${CMAKE_CURRENT_BINARY_DIR}/cmake/PETAce-DuetConfig.cmake)
set(DUET_TARGETS_FILENAME ${CMAKE_CURRENT_BINARY_DIR}/cmake/PETAce-DuetTargets.cmake)
set(DUET_CONFIG_VERSION_FILENAME ${CMAKE_CURRENT_BINARY_DIR}/cmake/PETAce-DuetConfigVersion.cmake)
set(DUET_CONFIG_H_FILENAME ${CMAKE_CURRENT_BINARY_DIR}/src/duet/util/config.h)
set(DUET_THIRDPARTY_DIR ${CMAKE_CURRENT_BINARY_DIR}/thirdparty)

# Installation tree
set(DUET_CONFIG_INSTALL_DIR ${CMAKE_INSTALL_LIBDIR}/cmake/PETAce-Duet-${DUET_VERSION_MAJOR}.${DUET_VERSION_MINOR})
set(DUET_INCLUDES_INSTALL_DIR ${CMAKE_INSTALL_INCLUDEDIR}/PETAce-${DUET_VERSION_MAJOR}.${DUET_VERSION_MINOR})
set(DUET_THIRDPARTY_INCLUDES_INSTALL_DIR ${DUET_INCLUDES_INSTALL_DIR}/thirdparty)

# Make the install target depend on the all target.
set(CMAKE_SKIP_INSTALL_ALL_DEPENDENCY OFF)

# Supported target operating systems are Linux and macOS.
if (NOT DEFINED LINUX)
    if (UNIX AND NOT APPLE AND NOT CYGWIN AND NOT MINGW)
        set(LINUX ON)
    endif()
endif()
if (UNIX AND APPLE)
    set(MACOS ON)
endif()
if (NOT LINUX AND NOT MACOS)
    message(FATAL_ERROR "Supported target operating systems are Linux and macOS")
endif()

# Only support x86_64 and arm64
set(CMAKE_REQUIRED_QUIET_OLD ${CMAKE_REQUIRED_QUIET})
set(CMAKE_REQUIRED_QUIET ON)
check_cxx_source_runs("
    #if defined(__aarch64__)
        int main() {
            return 0;
        }
    #else
        #error
    #endif
    "
    DUET_ARM64
)
check_cxx_source_runs("
    #if defined(__amd64)
        int main() {
            return 0;
        }
    #else
        #error
    #endif
    "
    DUET_AMD64
)
set(CMAKE_REQUIRED_QUIET ${CMAKE_REQUIRED_QUIET_OLD})
if (NOT DUET_AMD64 AND NOT DUET_ARM64)
    message(FATAL_ERROR "Supported target architectures are x86_64 and arm64")
endif()

add_compile_options(-msse4.2 -Wno-ignored-attributes -mavx)

set(DUET_ENABLE_GCOV_STR "Duet Enable gcov")
option(DUET_ENABLE_GCOV ${ENABLE_GCOV_STR} OFF)
message(STATUS "DUET_ENABLE_GCOV: ${DUET_ENABLE_GCOV_STR}")
if(CMAKE_BUILD_TYPE STREQUAL "Debug" AND DUET_ENABLE_GCOV)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
    set(CMAKE_CXX_LINK_EXECUTABLE "${CMAKE_CXX_LINK_EXECUTABLE} -fprofile-arcs -ftest-coverage -lgcov")
endif()

# CMAKE_CXX_LINK_EXECUTABLE
set(CMAKE_CXX_LINK_EXECUTABLE "${CMAKE_CXX_LINK_EXECUTABLE} -ldl -lrt")

#########################
# External dependencies #
#########################

# [OPTION] DUET_BUILD_DEPS (DEFAULT: ON)
# Download and build dependencies if set to ON.
# Look for dependencies using find_package, otherwise.
set(DUET_BUILD_DEPS_OPTION_STR "Automatically download and build unmet dependencies")
option(DUET_BUILD_DEPS ${DUET_BUILD_DEPS_OPTION_STR} ON)
message(STATUS "DUET_BUILD_DEPS: ${DUET_BUILD_DEPS}")

if(DUET_BUILD_DEPS)
    include(FetchContent)
    mark_as_advanced(FETCHCONTENT_BASE_DIR)
    mark_as_advanced(FETCHCONTENT_FULLY_DISCONNECTED)
    mark_as_advanced(FETCHCONTENT_UPDATES_DISCONNECTED)
    mark_as_advanced(FETCHCONTENT_QUIET)
endif()

# [OPTION] DUET_BUILD_SHARED_LIBS (DEFAULT: OFF)
# Build a shared library if set to ON.
set(DUET_BUILD_SHARED_LIBS_STR "Build shared library")
option(DUET_BUILD_SHARED_LIBS ${DUET_BUILD_SHARED_LIBS_STR} OFF)
message(STATUS "DUET_BUILD_SHARED_LIBS: ${DUET_BUILD_SHARED_LIBS}")

# Require Threads::Threads
if(NOT TARGET Threads::Threads)
    set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
    set(THREADS_PREFER_PTHREAD_FLAG TRUE)
    find_package(Threads REQUIRED)
    if(NOT Threads_FOUND)
        message(FATAL_ERROR "Threads: not found, please download and install manually")
    else()
        message(STATUS "Threads: found")
    endif()
endif()

# PETAce-Network::network
find_package(PETAce-Network 0.1 QUIET CONFIG)
if(PETAce-Network_FOUND)
    message(STATUS "PETAce-Network: found")
    if(PETAce-Network_STATIC_FOUND)
        set(network "PETAce-Network::network")
    else()
        set(network "PETAce-Network::network_shared")
    endif()
else()
    if(DUET_BUILD_DEPS)
        message(STATUS "PETAce-Network: download ...")
        duet_fetch_thirdparty_content(PETAce-Network)
        if(TARGET network)
            set(network "network")
        else()
            set(network "network_shared")
        endif()
        set(DUET_BUILD_NETWORK TRUE CACHE BOOL "" FORCE)
    else()
        message(FATAL_ERROR "PETAce-Network: not found, please download and install manually")
    endif()
endif()

# PETAce-Solo::solo
find_package(PETAce-Solo 0.1 QUIET CONFIG)
if(PETAce-Solo_FOUND)
    message(STATUS "PETAce-Solo: found")
    if(PETAce-Solo_STATIC_FOUND)
        set(solo "PETAce-Solo::solo")
    else()
        set(solo "PETAce-Solo::solo_shared")
    endif()
else()
    if(DUET_BUILD_DEPS)
        message(STATUS "PETAce-Solo: download ...")
        duet_fetch_thirdparty_content(PETAce-Solo)
        if(TARGET solo)
            set(solo "solo")
        else()
            set(solo "solo_shared")
        endif()
        set(DUET_BUILD_SOLO TRUE CACHE BOOL "" FORCE)
    else()
        message(FATAL_ERROR "PETAce-Solo: not found, please download and install manually")
    endif()
endif()

# PETAce-Verse::verse
find_package(PETAce-Verse 0.1 QUIET CONFIG)
if(PETAce-Verse_FOUND)
    message(STATUS "PETAce-Verse: found")
    if(PETAce-Verse_STATIC_FOUND)
        set(verse "PETAce-Verse::verse")
    else()
        set(verse "PETAce-Verse::verse_shared")
    endif()
else()
    if(DUET_BUILD_DEPS)
        message(STATUS "PETAce-Verse: download ...")
        duet_fetch_thirdparty_content(PETAce-Verse)
        if(TARGET verse)
            set(verse "verse")
        else()
            set(verse "verse_shared")
        endif()
        set(DUET_BUILD_VERSE TRUE CACHE BOOL "" FORCE)
    else()
        message(FATAL_ERROR "PETAce-Verse: not found, please download and install manually")
    endif()
endif()

# eigen
if(NOT TARGET Eigen3::Eigen)
    find_package(Eigen3 3.4.0 QUIET CONFIG)
    if(Eigen3_FOUND)
        message(STATUS "Eigen: found")
    else()
        if(DUET_BUILD_DEPS)
            message(STATUS "Eigen: download ...")
            duet_fetch_thirdparty_content(ExternalEigen)
            set(DUET_BUILD_EIGEN TRUE CACHE BOOL "" FORCE)
        else()
            message(FATAL_ERROR "Eigen: not found, please download and install manually")
        endif()
    endif()
endif()

####################
# DUET C++ library #
####################

# Add source files to library and header files to install
set(DUET_SOURCE_FILES "")
add_subdirectory(src/duet)

# Create the config file
configure_file(${DUET_CONFIG_H_IN_FILENAME} ${DUET_CONFIG_H_FILENAME})
install(
    FILES ${DUET_CONFIG_H_FILENAME}
    DESTINATION ${DUET_INCLUDES_INSTALL_DIR}/duet/util)

# Build only a static library
if(NOT DUET_BUILD_SHARED_LIBS)
    add_library(duet STATIC ${DUET_SOURCE_FILES})
    if(DUET_USE_CXX17)
        target_compile_features(duet PUBLIC cxx_std_17)
    else()
        target_compile_features(duet PUBLIC cxx_std_14)
    endif()
    target_include_directories(duet PUBLIC
        $<BUILD_INTERFACE:${DUET_INCLUDES_DIR}>
        $<INSTALL_INTERFACE:${DUET_INCLUDES_INSTALL_DIR}>)
    target_include_directories(duet PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/src/>)
    set_target_properties(duet PROPERTIES VERSION ${DUET_VERSION})
    set_target_properties(duet PROPERTIES OUTPUT_NAME petace_duet-${DUET_VERSION_MAJOR}.${DUET_VERSION_MINOR})

    if(DUET_BUILD_NETWORK)
        add_dependencies(duet ${network})
        target_include_directories(duet PUBLIC
            $<BUILD_INTERFACE:$<TARGET_PROPERTY:${network},INTERFACE_INCLUDE_DIRECTORIES>>
            $<INSTALL_INTERFACE:${DUET_INCLUDES_INSTALL_DIR}>)
        duet_combine_archives(duet ${network})
        set(DUET_CARRY_NETWORK TRUE)
    else()
        target_link_libraries(duet PUBLIC ${network})
        set(DUET_CARRY_NETWORK FALSE)
    endif()

    if(DUET_BUILD_SOLO)
        add_dependencies(duet ${solo})
        target_include_directories(duet PUBLIC
            $<BUILD_INTERFACE:$<TARGET_PROPERTY:${solo},INTERFACE_INCLUDE_DIRECTORIES>>
            $<INSTALL_INTERFACE:${DUET_INCLUDES_INSTALL_DIR}>)
    endif()
    target_link_libraries(duet PUBLIC ${solo})
    set(DUET_CARRY_SOLO FALSE)

    if(DUET_BUILD_VERSE)
        add_dependencies(duet ${verse})
        target_include_directories(duet PUBLIC
            $<BUILD_INTERFACE:$<TARGET_PROPERTY:${verse},INTERFACE_INCLUDE_DIRECTORIES>>
            $<INSTALL_INTERFACE:${DUET_INCLUDES_INSTALL_DIR}>)
    endif()
    target_link_libraries(duet PUBLIC ${verse})
    set(DUET_CARRY_VERSE FALSE)

    if(DUET_BUILD_EIGEN)
        add_dependencies(duet Eigen3::Eigen)
        target_include_directories(duet PUBLIC
            $<BUILD_INTERFACE:$<TARGET_PROPERTY:Eigen3::Eigen,INTERFACE_INCLUDE_DIRECTORIES>>
            $<INSTALL_INTERFACE:${DUET_THIRDPARTY_INCLUDES_INSTALL_DIR}>)
        set(DUET_CARRY_EIGEN TRUE)
    else()
        target_link_libraries(duet PUBLIC Eigen3::Eigen)
        set(DUET_CARRY_EIGEN FALSE)
    endif()

    target_link_libraries(duet PUBLIC Threads::Threads)

    install(TARGETS duet EXPORT PETAce-DuetTargets
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

# Build only a shared library
else()
    add_library(duet_shared SHARED ${DUET_SOURCE_FILES})
    if(DUET_USE_CXX17)
        target_compile_features(duet_shared PUBLIC cxx_std_17)
    else()
        target_compile_features(duet_shared PUBLIC cxx_std_14)
    endif()
    target_include_directories(duet_shared PUBLIC
        $<BUILD_INTERFACE:${DUET_INCLUDES_DIR}>
        $<INSTALL_INTERFACE:${DUET_INCLUDES_INSTALL_DIR}>)
    target_include_directories(duet_shared PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/src/>)
    set_target_properties(duet_shared PROPERTIES VERSION ${DUET_VERSION})
    set_target_properties(duet_shared PROPERTIES SOVERSION ${DUET_VERSION_MAJOR}.${DUET_VERSION_MINOR})
    set_target_properties(duet_shared PROPERTIES OUTPUT_NAME petace_duet)

    if(DUET_BUILD_NETWORK)
        target_include_directories(duet_shared PUBLIC
            $<BUILD_INTERFACE:$<TARGET_PROPERTY:${network},INTERFACE_INCLUDE_DIRECTORIES>>
            $<INSTALL_INTERFACE:${DUET_INCLUDES_INSTALL_DIR}>)
    endif()
    target_link_libraries(duet_shared PUBLIC ${network})
    set(DUET_CARRY_NETWORK FALSE)

    if(DUET_BUILD_SOLO)
        target_include_directories(duet_shared PUBLIC
            $<BUILD_INTERFACE:$<TARGET_PROPERTY:${solo},INTERFACE_INCLUDE_DIRECTORIES>>
            $<INSTALL_INTERFACE:${DUET_INCLUDES_INSTALL_DIR}>)
    endif()
    target_link_libraries(duet_shared PUBLIC ${solo})
    set(DUET_CARRY_SOLO FALSE)

    if(DUET_BUILD_VERSE)
        add_dependencies(duet_shared ${verse})
        target_include_directories(duet_shared PUBLIC
            $<BUILD_INTERFACE:$<TARGET_PROPERTY:${verse},INTERFACE_INCLUDE_DIRECTORIES>>
            $<INSTALL_INTERFACE:${DUET_INCLUDES_INSTALL_DIR}>)
    endif()
    target_link_libraries(duet_shared PUBLIC ${verse})
    set(DUET_CARRY_VERSE FALSE)

    if(DUET_BUILD_EIGEN)
        add_dependencies(duet_shared Eigen3::Eigen)
        target_include_directories(duet_shared PUBLIC
            $<BUILD_INTERFACE:$<TARGET_PROPERTY:Eigen3::Eigen,INTERFACE_INCLUDE_DIRECTORIES>>
            $<INSTALL_INTERFACE:${DUET_THIRDPARTY_INCLUDES_INSTALL_DIR}>)
        set(DUET_CARRY_EIGEN TRUE)
    else()
        target_link_libraries(duet_shared PUBLIC Eigen3::Eigen)
        set(DUET_CARRY_EIGEN FALSE)
    endif()

    target_link_libraries(duet_shared PUBLIC Threads::Threads)

    install(TARGETS duet_shared EXPORT PETAce-DuetTargets
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()

# Add standard alias targets for PETAce-Duet::duet and PETAce::duet_shared
if(TARGET duet)
    add_library(PETAce-Duet::duet ALIAS duet)
endif()
if(TARGET duet_shared)
    add_library(PETAce-Duet::duet_shared ALIAS duet_shared)
endif()

#################################
# Installation and CMake config #
#################################

# Create the CMake config file
include(CMakePackageConfigHelpers)
configure_package_config_file(
    ${DUET_CONFIG_IN_FILENAME} ${DUET_CONFIG_FILENAME}
    INSTALL_DESTINATION ${DUET_CONFIG_INSTALL_DIR})

# Install the export
install(
    EXPORT PETAce-DuetTargets
    NAMESPACE PETAce-Duet::
    DESTINATION ${DUET_CONFIG_INSTALL_DIR})

# Version file; we require exact version match for downstream
write_basic_package_version_file(
    ${DUET_CONFIG_VERSION_FILENAME}
    VERSION ${DUET_VERSION}
    COMPATIBILITY SameMinorVersion)

# Install config and module files
install(
    FILES
        ${DUET_CONFIG_FILENAME}
        ${DUET_CONFIG_VERSION_FILENAME}
    DESTINATION ${DUET_CONFIG_INSTALL_DIR})

# We export PETAce-DuetTargets from the build tree so it can be used by other projects
# without requiring an install.
export(
    EXPORT PETAce-DuetTargets
    NAMESPACE PETAce-Duet::
    FILE ${DUET_TARGETS_FILENAME})

# Install header files of dependencies if DUET_BUILD_DEPS is ON
if(DUET_BUILD_DEPS)
    # Insert dependencies here
    if(DUET_BUILD_NETWORK)
        install(CODE "execute_process(COMMAND ${CMAKE_COMMAND} --build ${network_BINARY_DIR} -t install)")
    endif()
    if(DUET_BUILD_SOLO)
        install(CODE "execute_process(COMMAND ${CMAKE_COMMAND} --build ${solo_BINARY_DIR} -t install)")
    endif()
    if(DUET_BUILD_VERSE)
        install(CODE "execute_process(COMMAND ${CMAKE_COMMAND} --build ${verse_BINARY_DIR} -t install)")
    endif()
    if(DUET_BUILD_EIGEN)
        install(
            DIRECTORY ${eigen_SOURCE_DIR}/Eigen
            DESTINATION ${DUET_THIRDPARTY_INCLUDES_INSTALL_DIR}
            COMPONENT Devel)
    endif()
endif()

####################
# DUET C++ example #
####################

# [option] DUET_BUILD_EXAMPLE
set(DUET_BUILD_EXAMPLE_OPTION_STR "Build C++ example for DUET")
option(DUET_BUILD_EXAMPLE ${DUET_BUILD_EXAMPLE_OPTION_STR} ON)
message(STATUS "DUET_BUILD_EXAMPLE: ${DUET_BUILD_EXAMPLE}")

if(DUET_BUILD_EXAMPLE)
    add_subdirectory(example)
endif()

#################
# DUET C++ test #
#################

# [option] DUET_BUILD_TEST
set(DUET_BUILD_TEST_OPTION_STR "Build C++ test for DUET")
option(DUET_BUILD_TEST ${DUET_BUILD_TEST_OPTION_STR} ON)
message(STATUS "DUET_BUILD_TEST: ${DUET_BUILD_TEST}")

if(DUET_BUILD_TEST)
    add_subdirectory(test)
    if(CMAKE_BUILD_TYPE STREQUAL "Debug" AND DUET_ENABLE_GCOV)
    add_custom_target(test_coverage
        COMMAND gcovr -r ${CMAKE_CURRENT_LIST_DIR} -f \"src\" -e \".+\(test\\.cpp\)\" --xml-pretty -o "${CMAKE_CURRENT_BINARY_DIR}/report/coverage.xml"
        WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
    endif()
endif()

##################
# DUET C++ bench #
##################

# [option] DUET_BUILD_BENCH
set(DUET_BUILD_BENCH_OPTION_STR "Build C++ benchmark for DUET")
option(DUET_BUILD_BENCH ${DUET_BUILD_BENCH_OPTION_STR} ON)
message(STATUS "DUET_BUILD_BENCH: ${DUET_BUILD_BENCH}")

if(DUET_BUILD_BENCH)
    add_subdirectory(bench)
endif()

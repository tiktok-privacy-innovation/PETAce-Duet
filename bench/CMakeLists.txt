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

project(DUETBench VERSION 0.3.0 LANGUAGES CXX)

# If not called from root CMakeLists.txt
if(NOT DEFINED DUET_BUILD_BENCH)
    set(DUET_BUILD_BENCH ON)

    find_package(PETAce-Duet 0.3.0 EXACT REQUIRED)

    add_compile_options(-msse4.2 -Wno-ignored-attributes -mavx)

    # Must define these variables and include macros
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
    set(DUET_THIRDPARTY_DIR ${CMAKE_CURRENT_BINARY_DIR}/thirdparty)
    set(THIRDPARTY_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/thirdparty)
    include(FetchContent)
    mark_as_advanced(FETCHCONTENT_BASE_DIR)
    mark_as_advanced(FETCHCONTENT_FULLY_DISCONNECTED)
    mark_as_advanced(FETCHCONTENT_UPDATES_DISCONNECTED)
    mark_as_advanced(FETCHCONTENT_QUIET)
    list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/../cmake)
    include(DuetCustomMacros)
else()
    set(THIRDPARTY_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/../thirdparty)
endif()

if(NOT DEFINED DUET_BUILD_DEPS)
    # [option] DUET_BUILD_DEPS (default: ON)
    # Download and build missing dependencies, throw error if disabled.
    set(DUET_BUILD_DEPS_OPTION_STR "Automatically download and build unmet dependencies")
    option(DUET_BUILD_DEPS ${DUET_BUILD_DEPS_OPTION_STR} ON)
endif()

# if DUET_BUILD_BENCH is ON, use GoogleBenchmark
if(DUET_BUILD_BENCH)

    if(DUET_BUILD_DEPS)
        # tclap
        message(STATUS "Tclap: download ...")
        duet_fetch_thirdparty_content(ExternalTclap)
    else()
        message(FATAL_ERROR "Tclap: not found, please download and install manually")
    endif()

    # glog::glog
    if(NOT TARGET glog::glog)
        find_package(glog 0.4.0 QUIET CONFIG)
        if(glog_FOUND)
            message(STATUS "glog: found")
        else()
            if(DUET_BUILD_DEPS)
                message(STATUS "glog: download ...")
                duet_fetch_thirdparty_content(ExternalGlog)
            else()
                message(FATAL_ERROR "glog: not found, please download and install manually")
            endif()
        endif()
    endif()

    # Add source files to bench
    set(DUET_BENCH_FILES
        ${CMAKE_CURRENT_LIST_DIR}/bench.cpp
        ${CMAKE_CURRENT_LIST_DIR}/duet_bench.cpp
    )

    set(CMAKE_CXX_LINK_EXECUTABLE "${CMAKE_CXX_LINK_EXECUTABLE} -ldl -lrt")
    add_executable(duet_bench ${DUET_BENCH_FILES})

    target_include_directories(duet_bench PRIVATE ${CMAKE_BINARY_DIR}/thirdparty/tclap-src/include)

    if(TARGET PETAce-Duet::duet)
        target_link_libraries(duet_bench PRIVATE PETAce-Duet::duet glog::glog)
    elseif(TARGET PETAce-Duet::duet_shared)
        target_link_libraries(duet_bench PRIVATE PETAce-Duet::duet_shared glog::glog)
    else()
        message(FATAL_ERROR "Cannot find target PETAce::duet or PETAce::duet_shared")
    endif()
endif()

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

FetchContent_Declare(
    seal
    GIT_REPOSITORY https://github.com/microsoft/SEAL.git
    GIT_TAG 206648d0e4634e5c61dcf9370676630268290b59
)

FetchContent_GetProperties(seal)

if(NOT seal_POPULATED)
    FetchContent_Populate(seal)

    set(SEAL_USE_CXX17 OFF CACHE BOOL "" FORCE)
    set(CMAKE_INSTALL_INCLUDEDIR_OLD ${CMAKE_INSTALL_INCLUDEDIR})
    set(CMAKE_INSTALL_INCLUDEDIR ${DUET_THIRDPARTY_INCLUDES_INSTALL_DIR})
    mark_as_advanced(FETCHCONTENT_SOURCE_DIR_SEAL)
    mark_as_advanced(FETCHCONTENT_UPDATES_DISCONNECTED_SEAL)

    add_subdirectory(
        ${seal_SOURCE_DIR}
        ${seal_BINARY_DIR}
        EXCLUDE_FROM_ALL)
    set(CMAKE_INSTALL_INCLUDEDIR ${CMAKE_INSTALL_INCLUDEDIR_OLD})
    unset(CMAKE_INSTALL_INCLUDEDIR_OLD)
endif()

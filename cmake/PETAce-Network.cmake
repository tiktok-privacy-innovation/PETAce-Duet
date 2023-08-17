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
    network
    GIT_REPOSITORY git@github.com:tiktok-privacy-innovation/PETAce-Network.git
    GIT_TAG        096245f85cec4d0065488225bf1b2f1dc979efb2 # v0.1.0
)
FetchContent_GetProperties(network)

if(NOT network_POPULATED)
    FetchContent_Populate(network)

    set(NETWORK_BUILD_SHARED_LIBS ${DUET_BUILD_SHARED_LIBS} CACHE BOOL "" FORCE)
    set(NETWORK_BUILD_TEST OFF CACHE BOOL "" FORCE)
    set(NETWORK_BUILD_EXAMPLE OFF CACHE BOOL "" FORCE)

    mark_as_advanced(FETCHCONTENT_SOURCE_DIR_NETWORK)
    mark_as_advanced(FETCHCONTENT_UPDATES_DISCONNECTED_NETWORK)

    add_subdirectory(
        ${network_SOURCE_DIR}
        ${network_SOURCE_DIR}/../network-build
        EXCLUDE_FROM_ALL)
endif()

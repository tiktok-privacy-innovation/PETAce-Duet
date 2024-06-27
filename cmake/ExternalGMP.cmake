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

include(ExternalProject)

set(GMP_SOURCES_DIR ${DUET_THIRDPARTY_DIR}/gmp)
set(GMP_INSTALL_DIR "${GMP_SOURCES_DIR}/install")
set(GMP_INCLUDE_DIR "${GMP_INSTALL_DIR}/include")
set(GMP_NAME "gmp")

include(ProcessorCount)
ProcessorCount(NUM_OF_PROCESSOR)

if((NOT DEFINED GMP_URL) OR (NOT DEFINED GMP_VER))
  message(STATUS "use pre defined download url")
  set(GMP_VER "gmp-6.3.0" CACHE STRING "" FORCE)
  set(GMP_URL "https://gmplib.org/download/gmp/${GMP_VER}.tar.bz2" CACHE STRING "" FORCE)
endif()

ExternalProject_Add(
  extern_gmp
  PREFIX            ${GMP_SOURCES_DIR}
  DOWNLOAD_COMMAND  wget --no-check-certificate ${GMP_URL} -c -q -O ${GMP_NAME}.tar.bz2
                    && tar -xvf ${GMP_NAME}.tar.bz2
  SOURCE_DIR      ${GMP_SOURCES_DIR}/src/${GMP_VER}
  CONFIGURE_COMMAND ./configure CFLAGS=-fPIC CXXFLAGS=-fPIC --prefix=${GMP_INSTALL_DIR}
                    --enable-cxx --with-pic
  BUILD_COMMAND     make -j ${NUM_OF_PROCESSOR}
  INSTALL_COMMAND   make install
  BUILD_IN_SOURCE   1
)

set(GMP_C_STATIC_LIBRARY "${GMP_INSTALL_DIR}/lib/libgmp.a")
set(GMP_CXX_STATIC_LIBRARY "${GMP_INSTALL_DIR}/lib/libgmpxx.a")

if(LINUX)
    set(GMP_C_SHARED_LIBRARY "${GMP_INSTALL_DIR}/lib/libgmp.so")
    set(GMP_CXX_SHARED_LIBRARY "${GMP_INSTALL_DIR}/lib/libgmpxx.so")
elseif(MACOS)
    set(GMP_C_SHARED_LIBRARY "${GMP_INSTALL_DIR}/lib/libgmp.dylib")
    set(GMP_CXX_SHARED_LIBRARY "${GMP_INSTALL_DIR}/lib/libgmpxx.dylib")
endif()

add_library(gmp STATIC IMPORTED GLOBAL)
set_property(TARGET gmp PROPERTY IMPORTED_LOCATION ${GMP_C_STATIC_LIBRARY})

add_library(gmp_shared SHARED IMPORTED GLOBAL)
set_property(TARGET gmp_shared PROPERTY IMPORTED_LOCATION ${GMP_C_SHARED_LIBRARY})

add_library(gmpxx STATIC IMPORTED GLOBAL)
set_property(TARGET gmpxx PROPERTY IMPORTED_LOCATION ${GMP_CXX_STATIC_LIBRARY})

add_library(gmpxx_shared SHARED IMPORTED GLOBAL)
set_property(TARGET gmpxx_shared PROPERTY IMPORTED_LOCATION ${GMP_CXX_SHARED_LIBRARY})

add_dependencies(gmp extern_gmp)
add_dependencies(gmpxx extern_gmp)
add_dependencies(gmp_shared extern_gmp)
add_dependencies(gmpxx_shared extern_gmp)

include_directories(${GMP_INCLUDE_DIR})

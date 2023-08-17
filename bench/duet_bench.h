// Copyright 2023 TikTok Pte. Ltd.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#pragma once

#include <cstdint>

#if (DUET_COMPILER == DUET_COMPILER_GCC)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#elif (DUET_COMPILER == DUET_COMPILER_CLANG)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wconversion"
#endif
#include "benchmark/benchmark.h"
#if (DUET_COMPILER == DUET_COMPILER_GCC)
#pragma GCC diagnostic pop
#elif (DUET_COMPILER == DUET_COMPILER_CLANG)
#pragma clang diagnostic pop
#endif

void bm_share_translation_name(benchmark::State& state, std::size_t rows, std::size_t cols);

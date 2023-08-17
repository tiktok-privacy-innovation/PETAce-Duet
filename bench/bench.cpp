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

#include "duet_bench.h"

#define PETACE_REG_BENCH(category, name, func, ...)                                                             \
    benchmark::RegisterBenchmark(                                                                               \
            (std::string(#category " / " #name)).c_str(), [=](benchmark::State& st) { func(st, __VA_ARGS__); }) \
            ->Unit(benchmark::kMicrosecond)                                                                     \
            ->Iterations(10);

int main(int argc, char** argv) {
    for (std::size_t i = 2; i < 8193; i *= 2) {
        for (std::size_t j = 1; j < 101; j *= 10) {
            PETACE_REG_BENCH(Shuffle, share_translation, bm_share_translation_name, i, j);
        }
    }

    benchmark::Initialize(&argc, argv);

    benchmark::RunSpecifiedBenchmarks();
}

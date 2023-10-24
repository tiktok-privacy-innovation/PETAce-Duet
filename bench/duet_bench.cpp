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

#include <array>
#include <cstdint>
#include <iostream>
#include <vector>

#include "duet/duet.h"
#include "duet/util/secret_shared_shuffle.h"

void bm_share_translation_name(benchmark::State& state, std::size_t rows, std::size_t cols) {
    std::cout << "rows:" << rows << " "
              << "cols:" << cols << std::endl;
    petace::duet::Permutation p(rows);
    petace::duet::ShareTranslation test(rows);
    for (auto _ : state) {
        std::vector<std::vector<petace::duet::OTChoice>> all_choice;
        test.active_phase_1(p, all_choice);
        std::vector<std::vector<std::vector<petace::duet::GGMTreeNode>>> all_levels_sums;
        petace::duet::Matrix<std::int64_t> a;
        petace::duet::Matrix<std::int64_t> b;
        test.passive_phase_1<std::int64_t>(cols, all_levels_sums, a, b);

        // mock_ot
        std::vector<std::vector<petace::duet::GGMTreeNode>> all_need_levels_sums(
                rows, std::vector<petace::duet::GGMTreeNode>(all_choice[0].size()));
        for (std::size_t i = 0; i < rows; ++i) {
            for (std::size_t j = 0; j < all_choice[0].size(); ++j) {
                if (all_choice[i][j] == 0) {
                    all_need_levels_sums[i][j] = all_levels_sums[i][j][0];
                } else {
                    all_need_levels_sums[i][j] = all_levels_sums[i][j][1];
                }
            }
        }

        petace::duet::Matrix<std::int64_t> delta;
        test.active_phase_2<std::int64_t>(cols, p, all_need_levels_sums, delta);
    }
}

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

#include "duet/st_generator/st_generator.h"

#include "duet/util/secret_shared_shuffle.h"

namespace petace {
namespace duet {

void STGenerator::get_st(std::size_t rows, std::size_t cols, const std::shared_ptr<network::Network>& net,
        Permutation& p, PlainMatrix<std::int64_t>& delta) {
    ShareTranslation st_generator(rows);
    std::vector<std::vector<OTChoice>> all_choice;
    st_generator.active_phase_1(p, all_choice);
    std::size_t log2_rows = all_choice[0].size();
    std::vector<OTChoice> ot_choice(rows * log2_rows);
    for (std::size_t i = 0; i < rows; ++i) {
        for (std::size_t j = 0; j < log2_rows; ++j) {
            ot_choice[i * log2_rows + j] = all_choice[i][j];
        }
    }
    std::vector<GGMTreeNode> ot_msgs;

    ot_->get_standard_ot<GGMTreeNode>(net, ot_choice, ot_msgs);

    std::vector<std::vector<GGMTreeNode>> all_need_levels_sums(rows, std::vector<GGMTreeNode>(log2_rows));
    for (std::size_t i = 0; i < rows; ++i) {
        for (std::size_t j = 0; j < log2_rows; ++j) {
            all_need_levels_sums[i][j] = ot_msgs[i * log2_rows + j];
        }
    }
    st_generator.active_phase_2<std::int64_t>(cols, p, all_need_levels_sums, delta);
}

void STGenerator::get_st(std::size_t rows, std::size_t cols, const std::shared_ptr<network::Network>& net,
        PlainMatrix<std::int64_t>& a, PlainMatrix<std::int64_t>& b) {
    ShareTranslation st_generator(rows);
    std::vector<std::vector<std::vector<GGMTreeNode>>> all_levels_sums;
    st_generator.passive_phase_1<std::int64_t>(cols, all_levels_sums, a, b);
    std::size_t log2_rows = all_levels_sums[0].size();
    std::vector<std::vector<GGMTreeNode>> ot_msgs(rows * log2_rows, std::vector<GGMTreeNode>(2));
    for (std::size_t i = 0; i < rows; ++i) {
        for (std::size_t j = 0; j < log2_rows; ++j) {
            ot_msgs[i * log2_rows + j] = all_levels_sums[i][j];
        }
    }
    ot_->get_standard_ot<GGMTreeNode>(net, ot_msgs);
}

}  // namespace duet
}  // namespace petace

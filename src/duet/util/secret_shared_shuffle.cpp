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

#include "duet/util/secret_shared_shuffle.h"

#include "duet/util/common.h"

namespace petace {
namespace duet {

// solo::AES DoublePrg::fixed_key_aes;
std::shared_ptr<solo::Hash> DoublePrg::fixed_hash = solo::Hash::create(solo::HashScheme::SHA_256);

void DoublePrg::gen(const PRGSeed& in, PRGSeed& first, PRGSeed& second) {
    // refer to https://eprint.iacr.org/2019/1340 & https://eprint.iacr.org/2019/074
    // G(x) = H(1 ^ x) || H(2 ^ x)
    // H(x) = fixed_aes(x) ^ x
    block one = in;
    block two = in;
    one[1] ^= 1;
    two[1] ^= 2;
    // TODO(@yindong):  change hash to aes.

    DoublePrg::fixed_hash->compute(
            reinterpret_cast<solo::Byte*>(&one), sizeof(block), reinterpret_cast<solo::Byte*>(&first), sizeof(block));
    DoublePrg::fixed_hash->compute(
            reinterpret_cast<solo::Byte*>(&two), sizeof(block), reinterpret_cast<solo::Byte*>(&second), sizeof(block));
    first ^= one;
    second ^= two;
}

ObliviousPuncturedVector::ObliviousPuncturedVector(std::size_t n) : n_(n) {
    log2n_ = ceil_log2(n);
}

std::size_t ObliviousPuncturedVector::get_ot_num() {
    return log2n_;
}

void ObliviousPuncturedVector::pos_to_choice(std::size_t pos, std::vector<GGMTreePuncChoice>& choice) {
    choice.resize(log2n_);
    for (std::size_t i = log2n_; i != 0; --i) {
        choice[i - 1] = pos % 2;
        pos = pos >> 1;
    }
}

void ObliviousPuncturedVector::active_phase_1(std::size_t pos, std::vector<OTChoice>& ot_choice) {
    if (pos >= n_) {
        throw std::invalid_argument("pos should less than n");
    }
    pos_to_choice(pos, ot_choice);
    for (std::size_t i = 0; i < log2n_; ++i) {
        ot_choice[i] = 1 ^ ot_choice[i];
    }
}

void ObliviousPuncturedVector::expand(std::size_t layer_now, std::vector<GGMTreeNode>& need_level) {
    std::size_t layer_items_num = 1;
    layer_items_num = layer_items_num << layer_now;
    GGMTreeNode first;
    GGMTreeNode second;
    for (std::size_t i = layer_items_num; i != 0; --i) {
        DoublePrg::gen(need_level[i - 1], first, second);
        need_level[2 * (i - 1)] = first;
        need_level[2 * (i - 1) + 1] = second;
    }
}

void ObliviousPuncturedVector::passive_phase_1(
        std::vector<std::vector<GGMTreeNode>>& levels_sums, std::vector<GGMTreeNode>& leaves) {
    levels_sums.resize(log2n_);
    leaves.resize(1ULL << log2n_);
    GGMTreeNode zero = {0, 0};
    // set root
    leaves[0] = read_block_from_dev_urandom();
    for (std::size_t i = 0; i < log2n_; ++i) {
        expand(i, leaves);
        levels_sums[i].resize(2);
        levels_sums[i][0] = zero;
        levels_sums[i][1] = zero;
        for (std::size_t j = 0; j < (1ULL << i); ++j) {
            levels_sums[i][0] ^= leaves[j * 2];
            levels_sums[i][1] ^= leaves[j * 2 + 1];
        }
    }
    leaves.erase(leaves.begin() + n_, leaves.end());
}

void ObliviousPuncturedVector::active_phase_2(
        std::size_t pos, const std::vector<GGMTreeNode>& need_levels_sums, std::vector<GGMTreeNode>& leaves) {
    if (pos >= n_) {
        throw std::invalid_argument("pos should less than n");
    }
    std::vector<GGMTreePuncChoice> choice;
    pos_to_choice(pos, choice);

    leaves.resize(1ULL << log2n_);
    GGMTreeNode zero = {0, 0};
    GGMTreeNode punc_tmp_node;
    // set root zero
    leaves[0] = zero;
    std::size_t index_tmp = 0;
    std::size_t choice_tmp;
    std::size_t level_sum_items_num;
    for (std::size_t i = 0; i < log2n_; ++i) {
        expand(i, leaves);
        level_sum_items_num = (1ULL << i);
        // set layer punc now is 0;
        leaves[2 * index_tmp] = zero;
        leaves[2 * index_tmp + 1] = zero;
        choice_tmp = 1 ? choice[i] == 0 : 0;
        punc_tmp_node = zero;
        for (std::size_t j = 0; j < level_sum_items_num; ++j) {
            punc_tmp_node ^= leaves[2 * j + choice_tmp];
        }
        punc_tmp_node ^= need_levels_sums[i];
        leaves[2 * index_tmp + choice_tmp] ^= punc_tmp_node;
        index_tmp = index_tmp * 2 + choice[i];
    }
    leaves.erase(leaves.begin() + n_, leaves.end());
}

ShareTranslation::ShareTranslation(std::size_t n) : n_(n) {
    opv_impl_ = std::make_unique<ObliviousPuncturedVector>(n);
}

void ShareTranslation::active_phase_1(const Permutation& p, std::vector<std::vector<OTChoice>>& all_ot_choice) {
    if (n_ != p.size()) {
        throw std::invalid_argument("n shuold equal permutation's size");
    }
    all_ot_choice.resize(n_);
    for (std::size_t i = 0; i < n_; ++i) {
        opv_impl_->active_phase_1(p[i], all_ot_choice[i]);
    }
}

void ShareTranslation::passive_phase_1(std::vector<std::vector<std::vector<GGMTreeNode>>>& all_levels_sums,
        std::vector<std::vector<GGMTreeNode>>& all_leaves_matrix) {
    all_levels_sums.resize(n_);
    all_leaves_matrix.resize(n_);
    for (std::size_t i = 0; i < n_; ++i) {
        opv_impl_->passive_phase_1(all_levels_sums[i], all_leaves_matrix[i]);
    }
}

void ShareTranslation::active_phase_2(const Permutation& p,
        const std::vector<std::vector<GGMTreeNode>>& all_need_levels_sums,
        std::vector<std::vector<GGMTreeNode>>& all_leaves_matrix) {
    all_leaves_matrix.resize(n_);
    for (std::size_t i = 0; i < n_; ++i) {
        opv_impl_->active_phase_2(p[i], all_need_levels_sums[i], all_leaves_matrix[i]);
    }
}

}  // namespace duet
}  // namespace petace

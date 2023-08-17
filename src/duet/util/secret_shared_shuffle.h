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

#include <memory>
#include <vector>

#include "solo/hash.h"
#include "solo/prng.h"

#include "duet/util/defines.h"
#include "duet/util/permutation.h"
#include "duet/util/prng.h"

namespace petace {
namespace duet {

/**
 * @brief Used in secret shared shuffle. For GGM tree.
 * @see Refer to https://eprint.iacr.org/2019/1340
 */
class DoublePrg {
public:
    // maybe global?
    static std::shared_ptr<solo::Hash> fixed_hash;

    /**
     * @brief PRG n -> 2n
     *
     * A PRG generate double length output of input. For GGM tree.
     *
     * @param[in] in input seed.
     * @param[out] first first part of output.
     * @param[out] second second part of output.
     */
    static void gen(const PRGSeed& in, PRGSeed& first, PRGSeed& second);
};

/**
 * @brief Implement of oblivious punctured vector.
 *
 * The interface is designed according to the way OT operates.
 * Active party as OT receiver, get punctured vector.
 * Passive party as OT sender, get whole vector.
 *
 * @see Refer to https://eprint.iacr.org/2019/1340
 */
class ObliviousPuncturedVector {
public:
    /**
     * @brief Constructor.
     *
     * Will be called in both party. n must same in two partys.
     *
     * @param[in] n Size of vector.
     */
    explicit ObliviousPuncturedVector(std::size_t n);

    ~ObliviousPuncturedVector() = default;

    std::size_t get_ot_num();

    /**
     * @brief Preparing for the ot choice corresponds to position.
     *
     * The inverse of OT choice also denote the path in the GGM tree from root to leaves.
     * For example, n = 12, pos is 9, OT choice is {0, 1, 1, 0};
     *
     * @param[in] pos Which position to punctured.
     * @param[out] ot_choice A 0/1 vector of length log2n, i.e. GGM tree's deepth.
     * @throws std::invalid_argument if pos not in [0, n)
     */
    void active_phase_1(std::size_t pos, std::vector<OTChoice>& ot_choice);

    /**
     * @brief Generate the whole vector and prepare OT messages.
     *
     * One OT message is the odd nodes sums and even nodes sums in the same level.
     * For example, in level 2 [0, 1, 2, 3] -> OT message: (tree[2][0] + tree[2][2],
     * tree[2][1] + tree[2][3])
     *
     * @param[out] levels_sums OT messages, i.e. sums of different levels of the GGM tree.
     * @param[out] leaves Leaves of GGM tree, i.e. the whole vector.
     */
    void passive_phase_1(std::vector<std::vector<GGMTreeNode>>& levels_sums, std::vector<GGMTreeNode>& leaves);

    /**
     * @brief Generate the punctured vector after OT.
     *
     * Using OT choice generated in active phase 1 and executing OT between both party
     * to get the need_levels_sums.
     *
     * @param[in] pos Which position to punctured.
     * @param[in] need_levels_sums OT results.
     * @param[out] leaves Leaves of punctured GGM tree, i.e. the punctured vector.
     * @throws std::invalid_argument if pos not in [0, n)
     */
    void active_phase_2(
            std::size_t pos, const std::vector<GGMTreeNode>& need_levels_sums, std::vector<GGMTreeNode>& leaves);

private:
    // noncopyable todo: maybe use boost
    ObliviousPuncturedVector(const ObliviousPuncturedVector&) = delete;

    void operator=(const ObliviousPuncturedVector&) = delete;

    // need malloc space before this function
    void expand(std::size_t layer_now, std::vector<GGMTreeNode>& need_level);

    void pos_to_choice(std::size_t pos, std::vector<GGMTreePuncChoice>& choice);

    std::size_t n_;
    std::size_t log2n_;
};

/**
 * @brief Implement of share translation.
 *
 * The interface is designed according to the way OT operates.
 * Active party as OT receiver, get punctured matrix.
 * Passive party as OT sender, get whole matrix.
 * The oblivious punctured vector protocol will be called n times.
 *
 * @see Refer to https://eprint.iacr.org/2019/1340
 */
class ShareTranslation {
public:
    /**
     * @brief Constructor.
     *
     * Will be called in both party. n must same in two partys.
     *
     * @param[in] n Size of permutation.
     */
    explicit ShareTranslation(std::size_t n);

    /**
     * @brief Preparing for the ot choice corresponds to permutation.
     *
     * Calling oblivious punctured vector protocol n times. Each use a position
     * in the permutation as pos for the oblivious punctured vector protocol.
     *
     * @param[in] p Permutation, usually is a ranom permutaion.
     * @param[out] all_ot_choice All ot choice generated by calling oblivious punctured vector protocol.
     */
    void active_phase_1(const Permutation& p, std::vector<std::vector<OTChoice>>& all_ot_choice);

    /**
     * @brief Generate the whole matrix and prepare OT messages.
     *
     * Calling oblivious punctured vector protocol n times.
     * Intermediate results.
     *
     * @param[out] all_levels_sums The OT messages needed by calling oblivious punctured vector protocol.
     * @param[out] all_leaves_matrix The whole matrix.
     */
    void passive_phase_1(std::vector<std::vector<std::vector<GGMTreeNode>>>& all_levels_sums,
            std::vector<std::vector<GGMTreeNode>>& all_leaves_matrix);

    /**
     * @brief Generate the whole matrix after OT.
     *
     * Using OT choice generated in active phase 1 and executing OT between both party
     * to get the need_levels_sums.
     * Intermediate results.
     *
     * @param[in] p Permutation, usually is a ranom permutaion.
     * @param[in] all_need_levels_sums OT results.
     * @param[out] all_leaves_matrix The whole matrix.
     */
    void active_phase_2(const Permutation& p, const std::vector<std::vector<GGMTreeNode>>& all_need_levels_sums,
            std::vector<std::vector<GGMTreeNode>>& all_leaves_matrix);

    /**
     * @brief Different from above function, this one Generate the final result.
     *
     * According to the protocol, we need usd a PRG extend the matrix and get the finally a&b.
     * p(a) - b = delta or p(a) xor b = delta.
     *
     * @param[in] cols The cols of the matrix(not the above matrix used in this protocol) needs to shuffle.
     * @param[in] is_boolean_share Defalut false. If it is ture, it generate p(a) xor b = delta.
     * @param[out] all_levels_sums The OT messages needed by calling oblivious punctured vector protocol.
     * @param[out] a Part of the tuple (p, a, b, delta)
     * @param[out] b Part of the tuple (p, a, b, delta)
     *
     */
    template <typename DataType>
    void passive_phase_1(std::size_t cols, std::vector<std::vector<std::vector<GGMTreeNode>>>& all_levels_sums,
            PlainMatrix<DataType>& a, PlainMatrix<DataType>& b, bool is_boolean_share = false) {
        std::vector<std::vector<GGMTreeNode>> all_leaves_matrix;
        passive_phase_1(all_levels_sums, all_leaves_matrix);

        std::vector<std::vector<std::vector<DataType>>> extend_leaves_matrix(n_);
        // get prng factory
        solo::PRNGFactory prng_factory(solo::PRNGScheme::AES_ECB_CTR);
        std::vector<solo::Byte> seed(sizeof(block));
        std::size_t cols_bytes_size = sizeof(DataType) * cols;

        for (std::size_t i = 0; i < n_; ++i) {
            extend_leaves_matrix[i].resize(n_);
            for (std::size_t j = 0; j < n_; ++j) {
                extend_leaves_matrix[i][j].resize(cols);
                memcpy(seed.data(), reinterpret_cast<solo::Byte*>(const_cast<GGMTreeNode*>(&all_leaves_matrix[i][j])),
                        sizeof(GGMTreeNode));

                auto prng = prng_factory.create(seed, cols_bytes_size);
                prng->generate(cols_bytes_size, reinterpret_cast<solo::Byte*>(extend_leaves_matrix[i][j].data()));
            }
        }
        std::vector<std::vector<DataType>> a_tmp(n_, std::vector<DataType>(cols, 0)),
                b_tmp(n_, std::vector<DataType>(cols, 0));
        for (std::size_t i = 0; i < n_; ++i) {
            for (std::size_t j = 0; j < n_; ++j) {
                for (std::size_t k = 0; k < cols; ++k) {
                    if (is_boolean_share) {
                        a_tmp[j][k] ^= extend_leaves_matrix[i][j][k];
                        b_tmp[i][k] ^= extend_leaves_matrix[i][j][k];
                    } else {
                        a_tmp[j][k] += extend_leaves_matrix[i][j][k];
                        b_tmp[i][k] += extend_leaves_matrix[i][j][k];
                    }
                }
            }
        }
        a.resize(n_, cols);
        b.resize(n_, cols);
        for (std::size_t i = 0; i < n_; ++i) {
            for (std::size_t j = 0; j < cols; ++j) {
                a(i * cols + j) = a_tmp[i][j];
                b(i * cols + j) = b_tmp[i][j];
            }
        }
    }

    /**
     * @brief Different from above function, this one Generate the final result.
     *
     * According to the protocol, we need usd a PRG extend the matrix and get the finally delta.
     * p(a) - b = delta or p(a) xor b = delta.
     *
     * @param[in] cols The cols of the matrix(not the above matrix used in this protocol) needs to shuffle.
     * @param[in] p Permutation, usually is a ranom permutaion.
     * @param[in] all_need_levels_sums OT results.
     * @param[in] is_boolean_share Defalut false. If it is ture, it generate p(a) xor b = delta.
     * @param[out] delta Part of the tuple (p, a, b, delta)
     *
     */
    template <typename DataType>
    void active_phase_2(std::size_t cols, const Permutation& p,
            const std::vector<std::vector<GGMTreeNode>>& all_need_levels_sums, PlainMatrix<DataType>& delta,
            bool is_boolean_share = false) {
        std::vector<std::vector<GGMTreeNode>> all_leaves_matrix;
        active_phase_2(p, all_need_levels_sums, all_leaves_matrix);

        std::vector<std::vector<std::vector<DataType>>> extend_leaves_matrix(n_);
        // get prng factory
        solo::PRNGFactory prng_factory(solo::PRNGScheme::AES_ECB_CTR);
        std::vector<solo::Byte> seed(sizeof(block));
        std::size_t cols_bytes_size = sizeof(DataType) * cols;

        for (std::size_t i = 0; i < n_; ++i) {
            extend_leaves_matrix[i].resize(n_);
            for (std::size_t j = 0; j < n_; ++j) {
                extend_leaves_matrix[i][j].resize(cols);
                memcpy(seed.data(), reinterpret_cast<solo::Byte*>(const_cast<GGMTreeNode*>(&all_leaves_matrix[i][j])),
                        sizeof(GGMTreeNode));

                auto prng = prng_factory.create(seed, cols_bytes_size);
                prng->generate(cols_bytes_size, reinterpret_cast<solo::Byte*>(extend_leaves_matrix[i][j].data()));
            }
        }

        std::vector<std::vector<DataType>> delta_tmp(n_, std::vector<DataType>(cols, 0));
        Permutation p_inv = p.inverse();
        for (std::size_t i = 0; i < n_; ++i) {
            for (std::size_t j = 0; j < n_; ++j) {
                for (std::size_t k = 0; k < cols; ++k) {
                    if (is_boolean_share) {
                        delta_tmp[i][k] ^= extend_leaves_matrix[i][j][k];
                        delta_tmp[p_inv[j]][k] ^= extend_leaves_matrix[i][j][k];
                    } else {
                        delta_tmp[i][k] -= extend_leaves_matrix[i][j][k];
                        delta_tmp[p_inv[j]][k] += extend_leaves_matrix[i][j][k];
                    }
                }
            }
        }
        delta.resize(n_, cols);
        for (std::size_t i = 0; i < n_; ++i) {
            for (std::size_t j = 0; j < cols; ++j) {
                delta(i * cols + j) = delta_tmp[i][j];
            }
        }
    }

private:
    // noncopyable todo: maybe use boost
    ShareTranslation(const ShareTranslation&) = delete;

    void operator=(const ShareTranslation&) = delete;

    std::size_t n_;
    std::unique_ptr<ObliviousPuncturedVector> opv_impl_ = nullptr;
};

/**
 * @brief Implement of Secret shared shuflle.
 *
 * Need to call twice if it is used to share.
 *
 * @see Refer to https://eprint.iacr.org/2019/1340
 */
class SecretSharedShuffle {
public:
    SecretSharedShuffle() = default;

    ~SecretSharedShuffle() = default;

    template <typename DataType>
    void passive_phase_1(const PlainMatrix<DataType>& x, const PlainMatrix<DataType>& a, PlainMatrix<DataType>& out) {
        out = x - a;
    }

    template <typename DataType>
    void active_phase_1(const Permutation& p, const PlainMatrix<DataType>& x_sub_a, const PlainMatrix<DataType>& delta,
            PlainMatrix<DataType>& out) {
        PlainMatrix<DataType> tmp = p.permute(x_sub_a);
        out = tmp + delta;
    }

private:
    // noncopyable todo: maybe use boost
    SecretSharedShuffle(const SecretSharedShuffle&) = delete;

    void operator=(const SecretSharedShuffle&) = delete;
};

}  // namespace duet
}  // namespace petace

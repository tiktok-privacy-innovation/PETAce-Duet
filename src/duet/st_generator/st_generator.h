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

#include "duet/ot_generator/ot_generator.h"
#include "duet/util/defines.h"
#include "duet/util/permutation.h"

namespace petace {
namespace duet {

/**
 * @brief Encapsulation of share translation.
 *
 * Used to duet to get share translation.
 *
 * @see Refer to https://eprint.iacr.org/2019/1340
 */
class STGenerator {
public:
    /**
     * @brief Constructor.
     *
     * @param[in] party Which party. 0/1.
     * @param[in] ot The OT generator.
     */
    STGenerator(std::size_t party, std::shared_ptr<OTGenerator> ot) : party_id_(party), ot_(ot) {
    }

    ~STGenerator() = default;

    /**
     * @brief Get p&delta.
     *
     * @param[in] rows Rows of the matrix need to permute.
     * @param[in] cols Cols of the matrix need to permute.
     * @param[in] net The network interface (e.g., PETAce-Network interface).
     * @param[in] p p of the tuple (p, delta, a, b);
     * @param[in] delta delta of the tuple (p, delta, a, b);
     */
    void get_st(std::size_t rows, std::size_t cols, const std::shared_ptr<network::Network>& net, Permutation& p,
            PlainMatrix<std::int64_t>& delta);

    /**
     * @brief Get a&b.
     *
     * @param[in] rows Rows of the matrix need to permute.
     * @param[in] cols Cols of the matrix need to permute.
     * @param[in] net The network interface (e.g., PETAce-Network interface).
     * @param[in] a a of the tuple (p, delta, a, b);
     * @param[in] b b of the tuple (p, delta, a, b);
     */
    void get_st(std::size_t rows, std::size_t cols, const std::shared_ptr<network::Network>& net,
            PlainMatrix<std::int64_t>& a, PlainMatrix<std::int64_t>& b);

private:
    STGenerator(const STGenerator&) = delete;

    void operator=(const STGenerator&) = delete;

    std::size_t party_id_ = 0;
    std::shared_ptr<OTGenerator> ot_ = nullptr;
};

}  // namespace duet
}  // namespace petace

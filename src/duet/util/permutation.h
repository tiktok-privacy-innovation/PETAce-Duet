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

#include <vector>

#include "duet/util/defines.h"

namespace petace {
namespace duet {

/**
 * @brief Permutation and related operations.
 *
 * The main operations are permute a vector like data, permute a permutation and
 * get the inverse of a permutaion.
 */
class Permutation {
public:
    Permutation() = delete;

    ~Permutation() = default;

    Permutation(Permutation&&) = default;

    Permutation& operator=(Permutation&&) = default;

    /**
     * @brief Generate the specified permutation.
     *
     * Data is like {1, 0, 2, 3}. Begin with 0 and less than n.
     *
     * @param[in] data Expression of permutation.
     * @throws std::invalid_argument if data is not invaild
     */
    explicit Permutation(const std::vector<std::size_t>& data);

    /**
     * @brief Generate a random permutation.
     *
     * @param[in] size The size of permutation.
     */
    explicit Permutation(std::size_t size);

    /**
     * @brief Generate a random permutation but it can set seed.
     *
     * @param[in] seed The random seed to generate random permutaion.
     * @param[in] size The size of permutation.
     */
    Permutation(const PRGSeed& seed, std::size_t size);

    const std::vector<std::size_t>& data() const;

    std::size_t size() const;

    std::size_t operator[](std::size_t index) const;

    /**
     * @brief Permute any vector like data.
     *
     * DataType must can be call with [].
     *
     * @param[in] in The input vector like data.
     * @return The result after permutation.
     */
    template <typename DataType>
    DataType permute(const DataType& in) const {
        DataType out;
        out.resize(in.size());
        for (std::size_t i = 0; i < size(); ++i) {
            out[i] = in[data_[i]];
        }
        return out;
    }

    /**
     * @brief Permute any vector like data with the inversed permutaion.
     *
     * DataType must can be call with [].
     *
     * @param[in] in The input vector like data.
     * @return The result after permutation with the inversed permutaion.
     */
    template <typename DataType>
    DataType inverse(const DataType& in) const {
        DataType out;
        out.resize(in.size());
        for (std::size_t i = 0; i < size(); ++i) {
            out[data_[i]] = in[i];
        }
        return out;
    }

    template <typename DataType>
    PlainMatrix<DataType> permute(const PlainMatrix<DataType>& in) const {
        PlainMatrix<DataType> out(in.rows(), in.cols());
        std::size_t rows = in.rows();
        for (std::size_t i = 0; i < rows; ++i) {
            out.row(i) = in.row(data_[i]);
        }
        return out;
    }

    /**
     * @brief Combine two permutaion.
     *
     * @param[in] in The permutation nedds to combine.
     * @return The combination of two permutation.
     */
    Permutation combine(const Permutation& other) const;

    ArithMatrix permute(const ArithMatrix& in) const;

    /**
     * @brief Get the inverse of the permutaion.
     *
     * @return The inverse of the permutaion.
     */
    Permutation inverse() const;

private:
    Permutation(const Permutation&) = delete;

    Permutation& operator=(const Permutation&) = delete;

    void gen_random_permutation(const PRGSeed& seed, std::size_t size);

    std::vector<std::size_t> data_{};
};

}  // namespace duet
}  // namespace petace

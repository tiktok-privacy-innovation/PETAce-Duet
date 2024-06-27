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

#include "solo/prng.h"

#include "duet/util/defines.h"

namespace petace {
namespace duet {

/**
 * @brief Used for duet protocol.
 *
 * Common means both party will get same random number.
 * Unique means both party will get independent random number.
 */
class PRNG {
public:
    PRNG() = delete;

    /**
     * @brief Constructor.
     *
     * Will be called in both party with same common_seed.
     *
     * @param[in] common_seed Seed of common PRNG which will generate same random number in both party.
     * @param[in] buff_size The size of the buffer. Default 256. Unit is 128 bits.
     */
    explicit PRNG(const block& common_seed, std::int64_t buff_size = 256);

    ~PRNG() = default;

    /**
     * Both party will get the same random number.
     *
     * @return A random std::int64_t number.
     */
    std::int64_t get_common_rand();

    /**
     * Each party will get different random number.
     *
     * @return A random int like number, depends on DataType.
     */
    template <typename DataType>
    DataType get_unique_rand() {
        if (unique_rand_idx + sizeof(DataType) > unique_rand_buff.size() * sizeof(block)) {
            refill_unique_buffer();
        }

        DataType ret = *(DataType*)((std::uint8_t*)unique_rand_buff.data() + unique_rand_idx);  // NOLINT
        unique_rand_idx += sizeof(DataType);
        return ret;
    }

    std::shared_ptr<solo::PRNG> get_unique_rand_gen() const {
        return unique_rand_gen;
    }

private:
    std::int64_t common_rand_idx = 0;
    std::int64_t unique_rand_idx = 0;
    std::shared_ptr<solo::PRNG> common_rand_gen = nullptr;
    std::shared_ptr<solo::PRNG> unique_rand_gen = nullptr;
    std::vector<block> common_rand_buff{};
    std::vector<block> unique_rand_buff{};

    void refill_common_buffer();

    void refill_unique_buffer();
};

}  // namespace duet
}  // namespace petace

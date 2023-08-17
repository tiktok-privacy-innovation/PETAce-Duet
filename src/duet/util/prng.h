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

#include <random>
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
    explicit PRNG(const block& common_seed, std::int64_t buff_size = 256) {
        solo::PRNGFactory prng_factory(solo::PRNGScheme::AES_ECB_CTR);

        // common_rand_gen_idx = 0;
        common_rand_buff.resize(buff_size);
        std::vector<solo::Byte> seed(sizeof(block));
        memcpy(seed.data(), reinterpret_cast<solo::Byte*>(const_cast<block*>(&common_seed)), sizeof(block));
        common_rand_gen = prng_factory.create(seed);

        // unique_rand_gen_idx = 0;
        unique_rand_buff.resize(buff_size);
        unique_rand_gen = prng_factory.create(sizeof(block));

        refill_common_buffer();
        refill_unique_buffer();
    }

    ~PRNG() = default;

    /**
     * Both party will get the same random number.
     *
     * @return A random std::int64_t number.
     */
    std::int64_t get_common_rand() {
        if (common_rand_idx + sizeof(std::int64_t) > common_rand_buff.size() * sizeof(block)) {
            refill_common_buffer();
        }
        std::int64_t ret = *(std::int64_t*)((std::uint8_t*)common_rand_buff.data() + common_rand_idx);  // NOLINT
        common_rand_idx += sizeof(std::int64_t);
        return ret;
    }

    /**
     * Each party will get different random number.
     *
     * @return A random std::int64_t number.
     */
    std::int64_t get_unique_rand() {
        if (unique_rand_idx + sizeof(std::int64_t) > unique_rand_buff.size() * sizeof(block)) {
            refill_unique_buffer();
        }

        std::int64_t ret = *(std::int64_t*)((std::uint8_t*)unique_rand_buff.data() + unique_rand_idx);  // NOLINT
        unique_rand_idx += sizeof(std::int64_t);
        return ret;
    }

    /**
     * Each party will get different random number.
     *
     * @return A random std::int8_t number.
     */
    std::int8_t get_unique_rand_int8() {
        if (unique_rand_idx + sizeof(std::int8_t) > unique_rand_buff.size() * sizeof(block)) {
            refill_unique_buffer();
        }

        std::int8_t ret = *(std::int8_t*)((std::uint8_t*)unique_rand_buff.data() + unique_rand_idx);  // NOLINT
        unique_rand_idx += sizeof(std::int8_t);
        return ret;
    }

private:
    std::int64_t common_rand_idx = 0;
    std::int64_t unique_rand_idx = 0;
    std::shared_ptr<solo::PRNG> common_rand_gen = nullptr;
    std::shared_ptr<solo::PRNG> unique_rand_gen = nullptr;
    std::vector<block> common_rand_buff{};
    std::vector<block> unique_rand_buff{};

    void refill_common_buffer() {
        common_rand_gen->generate(
                common_rand_buff.size() * sizeof(block), reinterpret_cast<solo::Byte*>(common_rand_buff.data()));
        common_rand_idx = 0;
    }

    void refill_unique_buffer() {
        unique_rand_gen->generate(
                unique_rand_buff.size() * sizeof(block), reinterpret_cast<solo::Byte*>(unique_rand_buff.data()));
        unique_rand_idx = 0;
    }
};

}  // namespace duet
}  // namespace petace

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

#include "duet/util/prng.h"

namespace petace {
namespace duet {

PRNG::PRNG(const block& common_seed, std::int64_t buff_size) {
    solo::PRNGFactory prng_factory(solo::PRNGScheme::AES_ECB_CTR);

    common_rand_buff.resize(buff_size);
    std::vector<solo::Byte> seed(sizeof(block));
    memcpy(seed.data(), reinterpret_cast<solo::Byte*>(const_cast<block*>(&common_seed)), sizeof(block));
    common_rand_gen = prng_factory.create(seed);

    unique_rand_buff.resize(buff_size);
    unique_rand_gen = prng_factory.create(sizeof(block));

    refill_common_buffer();
    refill_unique_buffer();
}

std::int64_t PRNG::get_common_rand() {
    if (common_rand_idx + sizeof(std::int64_t) > common_rand_buff.size() * sizeof(block)) {
        refill_common_buffer();
    }
    std::int64_t ret = *(std::int64_t*)((std::uint8_t*)common_rand_buff.data() + common_rand_idx);  // NOLINT
    common_rand_idx += sizeof(std::int64_t);
    return ret;
}

void PRNG::refill_common_buffer() {
    common_rand_gen->generate(
            common_rand_buff.size() * sizeof(block), reinterpret_cast<solo::Byte*>(common_rand_buff.data()));
    common_rand_idx = 0;
}

void PRNG::refill_unique_buffer() {
    unique_rand_gen->generate(
            unique_rand_buff.size() * sizeof(block), reinterpret_cast<solo::Byte*>(unique_rand_buff.data()));
    unique_rand_idx = 0;
}

}  // namespace duet
}  // namespace petace

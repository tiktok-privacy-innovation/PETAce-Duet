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

#include <cmath>

#include "solo/prng.h"

#include "solo/ahe_paillier.h"

#include "duet/util/consts.h"

namespace petace {
namespace duet {

inline block read_block_from_dev_urandom() {
    block ret;
    solo::PRNG::get_random_byte_array(sizeof(ret), reinterpret_cast<solo::Byte*>(&ret));
    return ret;
}

inline std::int64_t double_to_fixed(double input) {
    return static_cast<std::int64_t>(input * ((std::int64_t)1 << kFixedPointPrecision));
}

inline double fixed_to_double(std::int64_t input) {
    return static_cast<double>(input) / ((std::int64_t)1 << kFixedPointPrecision);
}

inline std::size_t ceil_log2(std::size_t in) {
    return static_cast<std::size_t>(std::ceil(std::log2(in)));
}

inline std::int64_t ipcl_bn_to_int64(const solo::ahepaillier::BigNum& in) {
    std::size_t length = in.DwordSize();
    if (length == 0) {
        return 0;
    }
    std::vector<std::uint32_t> data;
    in.num2vec(data);
    std::uint64_t value = data[0];
    if (length > 1) {
        value += static_cast<std::uint64_t>(data[1]) << 32;
    }
    return static_cast<std::int64_t>(value);
}

}  // namespace duet
}  // namespace petace

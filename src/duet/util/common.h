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

namespace petace {
namespace duet {

inline block read_block_from_dev_urandom() {
    block ret;
    solo::PRNG::get_random_byte_array(sizeof(ret), reinterpret_cast<solo::Byte*>(&ret));
    return ret;
}

inline std::int64_t float_to_fixed(double input) {
    return static_cast<std::int64_t>(input * ((std::int64_t)1 << 16));
}

inline double fixed_to_float(std::int64_t input) {
    return static_cast<double>(input) / ((std::int64_t)1 << 16);
}

inline std::size_t ceil_log2(std::size_t in) {
    return static_cast<std::size_t>(std::ceil(std::log2(in)));
}

}  // namespace duet
}  // namespace petace

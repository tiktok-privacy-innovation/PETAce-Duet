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

#include <iomanip>

#include "duet/util/defines.h"

namespace petace {
namespace duet {

template <typename DataType>
void print_vector(const std::vector<DataType>& in) {
    std::size_t width = (sizeof(DataType) * 8 / 10 + 1) * 3;
    for (std::size_t i = 0; i < in.size(); ++i) {
        std::cout << std::setiosflags(std::ios::left) << std::setw(width);
        if (sizeof(DataType) < sizeof(int)) {
            std::cout << static_cast<int>(in[i]);
        } else {
            std::cout << in[i];
        }
        if ((i + 1) % 8 == 0) {
            std::cout << std::endl;
        }
    }
    std::cout << std::endl;
}

inline void print_block(block& in) {
    for (std::size_t i = 0; i < 2; ++i) {
        std::cout << std::setiosflags(std::ios::left) << std::setw(21) << in[i];
    }
    std::cout << std::endl;
}

}  // namespace duet
}  // namespace petace

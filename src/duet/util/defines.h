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

#include <emmintrin.h>

#include <vector>

#include "Eigen/Dense"

#include "solo/ahe_paillier.h"

namespace petace {
namespace duet {

// some declarations
using block = __m128i;
using PRGSeed = block;
using OTLabel = block;
using OTChoice = std::int8_t;
using GGMTreeNode = block;
using GGMTreePuncChoice = std::int8_t;
using Buffer = std::vector<std::int8_t>;
using RegisterAddress = std::uint64_t;
using PublicDouble = double;
using PublicIndex = std::size_t;
using ByteVector = std::vector<solo::Byte>;

template <typename T>
using Matrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

}  // namespace duet
}  // namespace petace

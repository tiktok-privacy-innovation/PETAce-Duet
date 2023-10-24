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

#include <cstddef>

namespace petace {
namespace duet {

// parameters for cheetah secure comparison
// note: when setting up these parameters, we need to guarantee that the number of the blocks is a power of 2.
const std::size_t kBlockBitLength = 1;
const std::size_t kOTSize = std::size_t(1) << kBlockBitLength;
// arith share
const std::size_t kFixedPointPrecision = 16;
const std::size_t kPowDepth = 6;
// bool share
const std::size_t kKoggeStonePpaDepth = 6;
// ot
const std::size_t kDefaultBaseOtSizes = 128;
const std::size_t kDefaultExtOtSizes = 1024;
// triplet
const std::size_t kDefaulBooleanTripletBufferSize = 1024;
const std::size_t kDefaulArithmeticTripletBufferSize = 1024;
// parameters for HE
const std::size_t kPaillierKeySize = 2048;
const std::size_t kCipherByteSize = (kPaillierKeySize * 2 + 7) / 8;
const std::size_t kStatisticalLambda = 40;

}  // namespace duet
}  // namespace petace

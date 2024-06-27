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
#include <string>

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
const std::size_t kDefaultExtOtSizes = 8192;
// triple
const std::size_t kDefaulBooleanTripleBufferSize = 1024;
const std::size_t kDefaulArithmeticTripleBufferSize = 8192;
// parameters for FHE
const std::size_t kFHEBatchSize = 8192;
const std::size_t kPolyModulusDegree = 8192;
const std::size_t kCRTPrimeCount = 4;
const std::size_t kFHERandomBitLength = 168;
const std::string kTwoPowerSixtyFour = "18446744073709551616";
// parameters for HE
const std::size_t kPaillierKeySize = 2048;
const std::size_t kPaillierCipherSize = ((kPaillierKeySize + 7) / 8) * 2;
const std::size_t kStatisticalLambda = 40;
const std::size_t kPaillierThreads = 1;

// sigmoid
const double kSigmoidParmas[3] = {-0.018715, 0.24955, 0.4999};
// div
const double kTwoPointNine = 2.9142;
// max
const int64_t kMaxValue = 0x000000ffffffffff;
const int64_t kMinValue = 0x8fffffffffffffff;

}  // namespace duet
}  // namespace petace

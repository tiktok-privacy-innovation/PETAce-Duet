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

#include "network/network.h"
#include "solo/prng.h"

#include "duet/oblivious_transfer/oblivious_transfer.h"
#include "duet/util/prng.h"

namespace petace {
namespace duet {

/**
 * @brief Implementation of Beaver's arithmetic triplet generation protocol.
 *
 * This class implements the arithmetic triplet generation protocol based on OT.
 */
class ArithmeticTriple {
public:
    /**
     * @brief Constructor of the class.
     */
    ArithmeticTriple() = default;

    ~ArithmeticTriple() = default;

    /**
     * @brief Returns a Beaver Triples.
     *
     * @param[in] net The network instance (e.g., from PETAce-Network)
     * @return A random Beaver triplet
     */
    virtual std::vector<std::int64_t> get_rand_triple(const std::shared_ptr<network::Network>& net) = 0;

private:
    virtual void refill_rand_triple_buffer(const std::shared_ptr<network::Network>& net) = 0;

    virtual void gen_rand_triple(const std::shared_ptr<network::Network>& net) = 0;
};

}  // namespace duet
}  // namespace petace

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

#include "duet/ot_generator/ot_generator.h"

namespace petace {
namespace duet {

/**
 * @brief Implementation of Beaver's boolean triplet generation protocol.
 *
 * This class implements the boolean triplet generation protocol based on OT.
 */
class BooleanTriplet {
public:
    /**
     * @brief Constructor of the class.
     *
     * @param[in] net The network instance (e.g., from PETAce-Network)
     * @param[in] party The party id
     * @param[in] ot The OT instance
     */
    BooleanTriplet(
            const std::shared_ptr<network::Network>& net, std::size_t party_id, const std::shared_ptr<OTGenerator>& ot);

    ~BooleanTriplet() = default;

    /**
     * @brief Returns a Beaver triplets.
     *
     * @param[in] net The network instance (e.g., from PETAce-Network)
     * @param[in] party The party id
     * @return A random Beaver triplet
     */
    std::vector<std::int64_t> get_rand_triplet(const std::shared_ptr<network::Network>& net, std::size_t party);

private:
    void refill_rand_triplet_buffer(const std::shared_ptr<network::Network>& net, std::size_t party);

    void gen_rand_triplet(const std::shared_ptr<network::Network>& net, std::size_t party);

    int lsb_to_int64(const std::vector<std::int64_t>& in, std::vector<std::int64_t>& out);

    std::size_t rand_triplet_idx = 0;
    std::vector<std::vector<std::int64_t>> rand_triplet_buff{};
    std::size_t party_id_ = 0;
    std::shared_ptr<OTGenerator> ot_ = nullptr;
};

}  // namespace duet
}  // namespace petace

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

#include "duet/ot_generator/ot_generator.h"
#include "duet/util/prng.h"

namespace petace {
namespace duet {

/**
 * @brief Implementation of Beaver's arithmetic triplet generation protocol.
 *
 * This class implements the arithmetic triplet generation protocol based on OT.
 */
class ArithmeticTriplet {
public:
    ArithmeticTriplet() = default;

    ~ArithmeticTriplet() = default;

    /**
     * @brief Initializer of the class.
     *
     * @param[in] net The network instance (e.g., from PETAce-Network)
     * @param[in] party The party id
     * @param[in] prng The PRNG instance
     * @param[in] ot The OT instance
     */
    void initialize(const std::shared_ptr<network::Network>& net, std::size_t party, std::shared_ptr<PRNG> prng,
            const std::shared_ptr<OTGenerator>& ot) {
        party_id_ = party;
        ot_ = ot;
        prng_ = prng;
        rand_triplet_buff.resize(8192);
        refill_rand_triplet_buffer(net);
    }

    /**
     * @brief Returns a Beaver triplets.
     *
     * @param[in] net The network instance (e.g., from PETAce-Network)
     * @return A random Beaver triplet
     */
    std::vector<std::int64_t> get_rand_triplet(const std::shared_ptr<network::Network>& net) {
        if (rand_triplet_idx >= rand_triplet_buff.size()) {
            refill_rand_triplet_buffer(net);
        }
        std::vector<std::int64_t> ret;
        ret.emplace_back(rand_triplet_buff[rand_triplet_idx][0]);
        ret.emplace_back(rand_triplet_buff[rand_triplet_idx][1]);
        ret.emplace_back(rand_triplet_buff[rand_triplet_idx][2]);
        rand_triplet_idx++;
        return ret;
    }

private:
    void refill_rand_triplet_buffer(const std::shared_ptr<network::Network>& net) {
        gen_rand_triplet(net);
        rand_triplet_idx = 0;
        return;
    }

    void gen_rand_triplet(const std::shared_ptr<network::Network>& net) {
        std::size_t rtp_size = rand_triplet_buff.size();

        auto mul_share = [&](std::vector<std::int64_t>& a, std::vector<std::int64_t>& b, std::vector<std::int64_t>& q,
                                 std::vector<std::int64_t>& t) {
            a.resize(rtp_size);
            for (std::size_t i = 0; i < rtp_size; i++) {
                a[i] = prng_->get_unique_rand();
            }

            b.resize(rtp_size);
            q.resize(rtp_size);
            t.resize(rtp_size);

            std::vector<std::vector<std::int64_t>> send_msgs;
            std::vector<std::int64_t> recv_msgs;
            for (std::size_t i = 0; i < rtp_size; i++) {
                if (party_id_ == 0) {
                    ot_->get_correlated_ot(net, 64, a[i], send_msgs);
                    t[i] = 0;
                    for (std::size_t j = 0; j < 64; j++) {
                        t[i] += -(send_msgs[j][0] << j);
                    }
                } else {
                    std::vector<std::int8_t> choices;
                    ot_->get_correlated_ot(net, 64, choices, recv_msgs);
                    q[i] = 0;
                    for (std::size_t j = 0; j < 64; j++) {
                        q[i] += (recv_msgs[j] << j);
                    }
                    b[i] = 0;
                    for (std::size_t j = 0; j < 64; j++) {
                        b[i] |= static_cast<std::int64_t>(choices[j]) << j;
                    }
                }
            }
            return;
        };

        std::vector<std::int64_t> a0;
        std::vector<std::int64_t> b1;
        std::vector<std::int64_t> q1;
        std::vector<std::int64_t> t0;
        mul_share(a0, b1, q1, t0);

        std::vector<std::int64_t> b0;
        std::vector<std::int64_t> a1;
        std::vector<std::int64_t> t1;
        std::vector<std::int64_t> q0;
        mul_share(b0, a1, t1, q0);

        if (party_id_ == 0) {
            for (std::size_t i = 0; i < rtp_size; i++) {
                rand_triplet_buff[i].emplace_back(a0[i]);
                rand_triplet_buff[i].emplace_back(b0[i]);
                rand_triplet_buff[i].emplace_back(a0[i] * b0[i] + t0[i] + q0[i]);
            }
        } else {
            for (std::size_t i = 0; i < rtp_size; i++) {
                rand_triplet_buff[i].emplace_back(a1[i]);
                rand_triplet_buff[i].emplace_back(b1[i]);
                rand_triplet_buff[i].emplace_back(a1[i] * b1[i] + t1[i] + q1[i]);
            }
        }
        return;
    }

    std::size_t rand_triplet_idx = 0;
    std::vector<std::vector<std::int64_t>> rand_triplet_buff{};
    std::shared_ptr<OTGenerator> ot_ = nullptr;
    std::size_t party_id_ = 0;
    std::shared_ptr<PRNG> prng_ = nullptr;
};

}  // namespace duet
}  // namespace petace

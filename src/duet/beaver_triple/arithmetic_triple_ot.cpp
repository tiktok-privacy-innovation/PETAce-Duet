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

#include "duet/beaver_triple/arithmetic_triple_ot.h"

#include "duet/util/consts.h"

namespace petace {
namespace duet {

ArithmeticTripleOT::ArithmeticTripleOT(const std::shared_ptr<network::Network>& net, std::size_t party_id,
        std::shared_ptr<PRNG> prng, const std::shared_ptr<ObliviousTransfer>& ot)
        : party_id_(party_id), ot_(ot), prng_(prng) {
    rand_triple_buff_.resize(kDefaulArithmeticTripleBufferSize);
    refill_rand_triple_buffer(net);
}

std::vector<std::int64_t> ArithmeticTripleOT::get_rand_triple(const std::shared_ptr<network::Network>& net) {
    if (rand_triple_idx_ >= rand_triple_buff_.size()) {
        refill_rand_triple_buffer(net);
    }
    std::vector<std::int64_t> ret;
    ret.emplace_back(rand_triple_buff_[rand_triple_idx_][0]);
    ret.emplace_back(rand_triple_buff_[rand_triple_idx_][1]);
    ret.emplace_back(rand_triple_buff_[rand_triple_idx_][2]);
    rand_triple_idx_++;
    return ret;
}

void ArithmeticTripleOT::refill_rand_triple_buffer(const std::shared_ptr<network::Network>& net) {
    gen_rand_triple(net);
    rand_triple_idx_ = 0;
    return;
}

void ArithmeticTripleOT::gen_rand_triple(const std::shared_ptr<network::Network>& net) {
    std::size_t rtp_size = rand_triple_buff_.size();
    auto mul_share = [&](std::vector<std::int64_t>& a, std::vector<std::int64_t>& b, std::vector<std::int64_t>& q,
                             std::vector<std::int64_t>& t) {
        a.resize(rtp_size);
        for (std::size_t i = 0; i < rtp_size; i++) {
            a[i] = prng_->get_unique_rand<std::int64_t>();
        }

        std::vector<std::int64_t> delta(rtp_size * 64);
        for (std::size_t i = 0; i < rtp_size; i++) {
            for (std::size_t j = 0; j < 64; j++) {
                delta[i * 64 + j] = a[i];
            }
        }

        b.resize(rtp_size);
        q.resize(rtp_size);
        t.resize(rtp_size);

        std::vector<std::vector<std::int64_t>> send_msgs;
        std::vector<std::int64_t> recv_msgs;
        std::vector<std::int8_t> choices;
        if (party_id_ == 0) {
            ot_->get_batch_cot(net, 64 * rtp_size, delta, send_msgs);
        } else {
            ot_->get_batch_cot(net, 64 * rtp_size, choices, recv_msgs);
        }
        for (std::size_t i = 0; i < rtp_size; i++) {
            if (party_id_ == 0) {
                t[i] = 0;
                for (std::size_t j = 0; j < 64; j++) {
                    t[i] += -(send_msgs[i * 64 + j][0] << j);
                }
            } else {
                q[i] = 0;
                for (std::size_t j = 0; j < 64; j++) {
                    q[i] += (recv_msgs[i * 64 + j] << j);
                }
                b[i] = 0;
                for (std::size_t j = 0; j < 64; j++) {
                    b[i] |= static_cast<std::int64_t>(choices[i * 64 + j]) << j;
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
            rand_triple_buff_[i].emplace_back(a0[i]);
            rand_triple_buff_[i].emplace_back(b0[i]);
            rand_triple_buff_[i].emplace_back(a0[i] * b0[i] + t0[i] + q0[i]);
        }
    } else {
        for (std::size_t i = 0; i < rtp_size; i++) {
            rand_triple_buff_[i].emplace_back(a1[i]);
            rand_triple_buff_[i].emplace_back(b1[i]);
            rand_triple_buff_[i].emplace_back(a1[i] * b1[i] + t1[i] + q1[i]);
        }
    }
    return;
}

}  // namespace duet
}  // namespace petace

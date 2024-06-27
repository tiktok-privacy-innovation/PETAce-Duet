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

#include "duet/beaver_triple/boolean_triple.h"

#include "duet/util/consts.h"

namespace petace {
namespace duet {

BooleanTriple::BooleanTriple(const std::shared_ptr<network::Network>& net, std::size_t party_id,
        const std::shared_ptr<ObliviousTransfer>& ot)
        : party_id_(party_id), ot_(ot) {
    rand_triple_buff_.resize(kDefaulBooleanTripleBufferSize);
    refill_rand_triple_buffer(net, party_id_);
}

std::vector<std::int64_t> BooleanTriple::get_rand_triple(
        const std::shared_ptr<network::Network>& net, std::size_t party) {
    if (rand_triple_idx_ >= rand_triple_buff_.size()) {
        refill_rand_triple_buffer(net, party);
    }
    std::vector<std::int64_t> ret(3);
    ret[0] = rand_triple_buff_[rand_triple_idx_][0];
    ret[1] = rand_triple_buff_[rand_triple_idx_][1];
    ret[2] = rand_triple_buff_[rand_triple_idx_][2];
    rand_triple_idx_++;
    return ret;
}

void BooleanTriple::refill_rand_triple_buffer(const std::shared_ptr<network::Network>& net, std::size_t party) {
    gen_rand_triple(net, party);
    rand_triple_idx_ = 0;
    return;
}

void BooleanTriple::gen_rand_triple(const std::shared_ptr<network::Network>& net, std::size_t party) {
    std::size_t tp_buff_size = rand_triple_buff_.size();
    std::vector<std::int64_t> msgs0(64 * tp_buff_size);
    std::vector<std::int64_t> msgs1(64 * tp_buff_size);
    std::vector<std::int64_t> x0(tp_buff_size);
    std::vector<std::int64_t> x1(tp_buff_size);
    std::vector<std::int64_t> xa(tp_buff_size);
    std::vector<std::int64_t> a0(tp_buff_size);
    std::vector<std::int64_t> b1(tp_buff_size);
    std::vector<std::int64_t> u0(tp_buff_size);
    std::vector<std::int64_t> v0(tp_buff_size);

    if (party == 0) {
        for (std::size_t i = 0; i < 64 * tp_buff_size; i++) {
            std::vector<std::int64_t> msg;
            ot_->get_random_ot(net, msg);
            msgs0[i] = msg[0];
            msgs1[i] = msg[1];
        }
        lsb_to_int64(msgs0, x0);
        lsb_to_int64(msgs1, x1);
    } else {
        for (std::size_t i = 0; i < 64 * tp_buff_size; i++) {
            std::int64_t msg;
            std::int8_t choice;
            ot_->get_random_ot(net, choice, msg);
            msgs0[i] = msg;
            msgs1[i] = static_cast<std::int64_t>(choice);
        }
        lsb_to_int64(msgs0, xa);
        lsb_to_int64(msgs1, a0);
    }

    if (party == 0) {
        for (std::size_t i = 0; i < tp_buff_size; i++) {
            b1[i] = x0[i] ^ x1[i];
            v0[i] = x0[i];
        }
    } else {
        for (std::size_t i = 0; i < tp_buff_size; i++) {
            u0[i] = xa[i];
        }
    }

    std::vector<std::int64_t> a1(tp_buff_size), b0(tp_buff_size), u1(tp_buff_size), v1(tp_buff_size);

    if (party == 0) {
        for (std::size_t i = 0; i < 64 * tp_buff_size; i++) {
            std::vector<std::int64_t> msg;
            ot_->get_random_ot(net, msg);
            msgs0[i] = msg[0];
            msgs1[i] = msg[1];
        }
        lsb_to_int64(msgs0, x0);
        lsb_to_int64(msgs1, x1);
    } else {
        for (std::size_t i = 0; i < 64 * tp_buff_size; i++) {
            std::int64_t msg;
            std::int8_t choice;
            ot_->get_random_ot(net, choice, msg);
            msgs0[i] = msg;
            msgs1[i] = static_cast<std::int64_t>(choice);
        }
        lsb_to_int64(msgs0, xa);
        lsb_to_int64(msgs1, b0);
    }

    if (party == 0) {
        for (std::size_t i = 0; i < tp_buff_size; i++) {
            a1[i] = x0[i] ^ x1[i];
            v1[i] = x0[i];
        }
    } else {
        for (std::size_t i = 0; i < tp_buff_size; i++) {
            u1[i] = xa[i];
        }
    }

    for (std::size_t i = 0; i < tp_buff_size; i++) {
        rand_triple_buff_[i].resize(3);
        if (party == 0) {
            rand_triple_buff_[i][0] = b1[i];
            rand_triple_buff_[i][1] = a1[i];
            rand_triple_buff_[i][2] = (a1[i] & b1[i]) ^ v0[i] ^ v1[i];

        } else {
            rand_triple_buff_[i][0] = b0[i];
            rand_triple_buff_[i][1] = a0[i];
            rand_triple_buff_[i][2] = (a0[i] & b0[i]) ^ u0[i] ^ u1[i];
        }
    }
}

int BooleanTriple::lsb_to_int64(const std::vector<std::int64_t>& in, std::vector<std::int64_t>& out) {
    if ((in.size() % 64) != 0) {
        return -1;
    }

    std::size_t k = 0;
    for (std::size_t i = 0; i < in.size(); i += 64, k++) {
        out[k] = 0;
        for (std::size_t j = 0; j < 64; j++) {
            out[k] = (out[k] << 1) | (in[i + j] & 0x1);
        }
    }
    return 0;
}

}  // namespace duet
}  // namespace petace

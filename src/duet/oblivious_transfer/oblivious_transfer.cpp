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

#include "duet/oblivious_transfer/oblivious_transfer.h"

#include <stdexcept>

#include "solo/hash.h"

#include "duet/util/common.h"

namespace petace {
namespace duet {

ObliviousTransfer::ObliviousTransfer(
        const std::size_t party, const block seed, const std::size_t base_ot_sizes, const std::size_t ext_ot_sizes)
        : party_id_(party) {
    ex_choices_.resize(ext_ot_sizes / (sizeof(block) * 8));

    solo::PRNGFactory prng_factory(solo::PRNGScheme::AES_ECB_CTR);
    std::vector<solo::Byte> seed_byte(reinterpret_cast<solo::Byte*>(const_cast<block*>(&seed)),
            reinterpret_cast<solo::Byte*>(const_cast<block*>(&seed)) + sizeof(block));
    prng_ = prng_factory.create(seed_byte, sizeof(block));

    buf_sizes_ = ext_ot_sizes;
    verse::VerseParams params;
    params.base_ot_sizes = base_ot_sizes;
    params.ext_ot_sizes = ext_ot_sizes;

    if (party == 0) {
        base_ot_recver_ = verse::VerseFactory<petace::verse::BaseOtReceiver>::get_instance().build(
                verse::OTScheme::NaorPinkasReceiver, params);
        ex_ot_sender_ = verse::VerseFactory<petace::verse::OtExtSender>::get_instance().build(
                verse::OTScheme::IknpSender, params);
        base_ot_sender_ = verse::VerseFactory<petace::verse::BaseOtSender>::get_instance().build(
                verse::OTScheme::NaorPinkasSender, params);
        ex_ot_recver_ = verse::VerseFactory<petace::verse::OtExtReceiver>::get_instance().build(
                verse::OTScheme::IknpReceiver, params);
    } else {
        base_ot_sender_ = verse::VerseFactory<petace::verse::BaseOtSender>::get_instance().build(
                verse::OTScheme::NaorPinkasSender, params);
        ex_ot_recver_ = verse::VerseFactory<petace::verse::OtExtReceiver>::get_instance().build(
                verse::OTScheme::IknpReceiver, params);
        base_ot_recver_ = verse::VerseFactory<petace::verse::BaseOtReceiver>::get_instance().build(
                verse::OTScheme::NaorPinkasReceiver, params);
        ex_ot_sender_ = verse::VerseFactory<petace::verse::OtExtSender>::get_instance().build(
                verse::OTScheme::IknpSender, params);
    }
}

void ObliviousTransfer::initialize(const std::shared_ptr<network::Network>& net) {
    if (party_id_ == 0) {
        std::vector<block> base_recv_ots;
        block rand;
        prng_->generate(sizeof(block), reinterpret_cast<solo::Byte*>(&rand));
        base_choices_.emplace_back(rand);
        base_ot_recver_->receive(net, base_choices_, base_recv_ots);
        ex_ot_sender_->set_base_ots(base_choices_, base_recv_ots);

        std::vector<std::array<block, 2>> base_send_ots;
        base_ot_sender_->send(net, base_send_ots);
        ex_ot_recver_->set_base_ots(base_send_ots);
    } else {
        std::vector<std::array<block, 2>> base_send_ots;
        base_ot_sender_->send(net, base_send_ots);
        ex_ot_recver_->set_base_ots(base_send_ots);

        std::vector<block> base_recv_ots;
        block rand;
        prng_->generate(sizeof(block), reinterpret_cast<solo::Byte*>(&rand));
        base_choices_.emplace_back(rand);
        base_ot_recver_->receive(net, base_choices_, base_recv_ots);
        ex_ot_sender_->set_base_ots(base_choices_, base_recv_ots);
    }

    if (party_id_ == 0) {
        refill_buffer_sender(net);
        refill_buffer_receiver(net);
    } else {
        refill_buffer_receiver(net);
        refill_buffer_sender(net);
    }
}

// sender
void ObliviousTransfer::refill_buffer_sender(const std::shared_ptr<network::Network>& net) {
    ex_ot_sender_->send(net, send_msgs_);
    cur_idx_sender_ = 0;
}

// receiver
void ObliviousTransfer::refill_buffer_receiver(const std::shared_ptr<network::Network>& net) {
    for (std::size_t i = 0; i < ex_choices_.size(); i++) {
        prng_->generate(sizeof(block), reinterpret_cast<solo::Byte*>(&ex_choices_[i]));
    }
    ex_ot_recver_->receive(net, ex_choices_, recv_msgs_);
    cur_idx_recver_ = 0;
}

void ObliviousTransfer::get_correlated_ot(const std::shared_ptr<network::Network>& net, const std::size_t ot_size,
        const std::int64_t delta, std::vector<std::vector<std::int64_t>>& msgs) {
    msgs.resize(ot_size);
    for (std::size_t i = 0; i < ot_size; i++) {
        get_random_ot(net, msgs[i]);
    }

    std::vector<std::int64_t> y(ot_size);
    for (std::size_t i = 0; i < ot_size; i++) {
        y[i] = msgs[i][0] + msgs[i][1] + delta;
    }

    for (std::size_t i = 0; i < ot_size; i++) {
        msgs[i][1] = msgs[i][0] + delta;
    }

    net->send_data(y.data(), ot_size * sizeof(std::int64_t));
}

void ObliviousTransfer::get_correlated_ot(const std::shared_ptr<network::Network>& net, const std::size_t ot_size,
        std::vector<std::int8_t>& choices, std::vector<std::int64_t>& msgs) {
    choices.resize(ot_size);
    msgs.resize(ot_size);
    for (std::size_t i = 0; i < ot_size; i++) {
        get_random_ot(net, choices[i], msgs[i]);
    }

    std::vector<std::int64_t> y(ot_size);
    net->recv_data(y.data(), ot_size * sizeof(std::int64_t));

    for (std::size_t i = 0; i < ot_size; i++) {
        if (choices[i] != 0) {
            msgs[i] = y[i] - msgs[i];
        }
    }

    return;
}

void ObliviousTransfer::get_batch_cot(const std::shared_ptr<network::Network>& net, const std::size_t ot_size,
        const std::vector<std::int64_t>& delta, std::vector<std::vector<std::int64_t>>& msgs) {
    msgs.resize(ot_size);
    for (std::size_t i = 0; i < ot_size; i++) {
        get_random_ot(net, msgs[i]);
    }

    std::vector<std::int64_t> y(ot_size);
    for (std::size_t i = 0; i < ot_size; i++) {
        y[i] = msgs[i][0] + msgs[i][1] + delta[i];
    }

    for (std::size_t i = 0; i < ot_size; i++) {
        msgs[i][1] = msgs[i][0] + delta[i];
    }

    net->send_data(y.data(), ot_size * sizeof(std::int64_t));
}

void ObliviousTransfer::get_batch_cot(const std::shared_ptr<network::Network>& net, const std::size_t ot_size,
        std::vector<std::int8_t>& choices, std::vector<std::int64_t>& msgs) {
    choices.resize(ot_size);
    msgs.resize(ot_size);
    for (std::size_t i = 0; i < ot_size; i++) {
        get_random_ot(net, choices[i], msgs[i]);
    }

    std::vector<std::int64_t> y(ot_size);
    net->recv_data(y.data(), ot_size * sizeof(std::int64_t));

    for (std::size_t i = 0; i < ot_size; i++) {
        if (choices[i] != 0) {
            msgs[i] = y[i] - msgs[i];
        }
    }

    return;
}

// sender, only support 1-of-n
void ObliviousTransfer::get_nch1_ot(const std::shared_ptr<network::Network>& net, const std::size_t select_size,
        const std::vector<std::vector<std::int64_t>>& msgs) {
    std::size_t size = msgs.size();
    std::vector<std::int64_t> y(size * select_size);

    std::size_t log_select_size = ceil_log2(select_size);
    std::vector<std::vector<std::vector<block>>> ot_msgs(size);
    for (std::size_t i = 0; i < size; i++) {
        ot_msgs[i].resize(log_select_size);
        for (std::size_t j = 0; j < log_select_size; j++) {
            if (cur_idx_sender_ == buf_sizes_) {
                refill_buffer_sender(net);
                cur_idx_sender_ = 0;
            }
            ot_msgs[i][j].emplace_back(send_msgs_[cur_idx_sender_][0]);
            ot_msgs[i][j].emplace_back(send_msgs_[cur_idx_sender_][1]);
            cur_idx_sender_++;
        }
    }

    std::vector<std::int8_t> xor_choices(size * log_select_size);
    net->recv_data(xor_choices.data(), size * log_select_size);

    std::vector<std::vector<std::vector<block>>> xor_msgs(size);
    auto hash = solo::Hash::create(solo::HashScheme::SHA_256);
    for (std::size_t i = 0; i < size; i++) {
        xor_msgs[i].resize(log_select_size);
        for (std::size_t j = 0; j < log_select_size; j++) {
            xor_msgs[i][j].emplace_back(ot_msgs[i][j][xor_choices[i * log_select_size + j]]);
            xor_msgs[i][j].emplace_back(ot_msgs[i][j][xor_choices[i * log_select_size + j] ^ 1]);
        }
        for (std::size_t j = 0; j < select_size; j++) {
            y[i * select_size + j] = msgs[i][j];
            for (std::size_t k = 0; k < log_select_size; k++) {
                std::vector<solo::Byte> hash_in(reinterpret_cast<solo::Byte*>(&xor_msgs[i][k][(j >> k) & 0x1]),
                        reinterpret_cast<solo::Byte*>(&xor_msgs[i][k][(j >> k) & 0x1]) + sizeof(block));
                hash_in.emplace_back(static_cast<solo::Byte>(j));
                block hash_out;
                hash->compute(
                        hash_in.data(), sizeof(block) + 1, reinterpret_cast<solo::Byte*>(&hash_out), sizeof(block));
                y[i * select_size + j] ^= hash_out[0];
            }
        }
    }
    net->send_data(y.data(), size * select_size * sizeof(std::int64_t));

    return;
}

// receiver, only support 1-of-n
void ObliviousTransfer::get_nch1_ot(const std::shared_ptr<network::Network>& net, const std::size_t select_size,
        const std::vector<std::int8_t>& choices, std::vector<std::int64_t>& msgs) {
    std::size_t size = choices.size();
    std::vector<std::vector<block>> ot_msgs(size);

    std::size_t log_select_size = ceil_log2(select_size);
    std::vector<std::vector<std::int8_t>> rand_choices(size);
    for (std::size_t i = 0; i < size; i++) {
        ot_msgs[i].resize(log_select_size);
        rand_choices[i].resize(log_select_size);
        for (std::size_t j = 0; j < log_select_size; j++) {
            if (cur_idx_recver_ == buf_sizes_) {
                refill_buffer_receiver(net);
                cur_idx_recver_ = 0;
            }

            ot_msgs[i][j] = recv_msgs_[cur_idx_recver_];
            rand_choices[i][j] = static_cast<std::int8_t>(verse::bit_from_blocks(ex_choices_, cur_idx_recver_));
            cur_idx_recver_++;
        }
    }

    std::vector<std::int8_t> xor_choices(size * log_select_size);
    for (std::size_t i = 0; i < size; i++) {
        for (std::size_t j = 0; j < log_select_size; j++) {
            xor_choices[i * log_select_size + j] =
                    static_cast<std::int8_t>(rand_choices[i][j] ^ ((choices[i] >> j) & 0x1));
        }
    }
    net->send_data(xor_choices.data(), size * log_select_size);

    std::vector<std::int64_t> y(size * select_size);
    net->recv_data(y.data(), size * select_size * sizeof(std::int64_t));

    msgs.resize(size);
    auto hash = solo::Hash::create(solo::HashScheme::SHA_256);
    for (std::size_t i = 0; i < size; i++) {
        msgs[i] = y[i * select_size + choices[i]];
        for (std::size_t k = 0; k < log_select_size; k++) {
            std::vector<solo::Byte> hash_in(reinterpret_cast<solo::Byte*>(&ot_msgs[i][k]),
                    reinterpret_cast<solo::Byte*>(&ot_msgs[i][k]) + sizeof(block));
            hash_in.emplace_back(static_cast<solo::Byte>(choices[i]));
            block hash_out;
            hash->compute(hash_in.data(), sizeof(block) + 1, reinterpret_cast<solo::Byte*>(&hash_out), sizeof(block));
            msgs[i] ^= hash_out[0];
        }
    }

    return;
}

}  // namespace duet
}  // namespace petace

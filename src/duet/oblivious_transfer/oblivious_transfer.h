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

#include <array>
#include <memory>
#include <queue>
#include <string>
#include <vector>

#include "network/network.h"
#include "solo/prng.h"
#include "verse/base-ot/base_ot_receiver.h"
#include "verse/base-ot/base_ot_sender.h"
#include "verse/base-ot/naor-pinkas-ot/naor_pinkas_ot.h"
#include "verse/two-choose-one/iknp/iknp_ot_ext.h"
#include "verse/two-choose-one/ot_ext_receiver.h"
#include "verse/two-choose-one/ot_ext_sender.h"
#include "verse/util/common.h"
#include "verse/verse_factory.h"

#include "duet/util/consts.h"
#include "duet/util/defines.h"
#include "duet/util/io.h"

namespace petace {
namespace duet {

/**
 * @brief This class contains functions relevant to Oblivious transfer.
 *
 * It provide APIs to support random OT, 2-choose-1 OT and general n-choose-1 OT.
 */
class ObliviousTransfer {
public:
    ObliviousTransfer() = delete;

    /**
     * @brief Constructor of the class.
     *
     * @param[in] party The party id
     * @param[in] seed The seed for PRNG
     * @param[in] base_ot_sizes Optional paramater to initialize PETAce-Verse
     * @param[in] ext_ot_sizes Optional paramater to initialize PETAce-Verse
     */
    ObliviousTransfer(const std::size_t party, const block seed, const std::size_t base_ot_sizes = kDefaultBaseOtSizes,
            const std::size_t ext_ot_sizes = kDefaultExtOtSizes);

    ~ObliviousTransfer() = default;

    /**
     * @brief Initializer of the class.
     *
     * @param[in] net The network instance (e.g., from PETAce-Network)
     */
    void initialize(const std::shared_ptr<network::Network>& net);

    /**
     * @brief Refill the OT buffer for the sender
     *
     * @param[in] net The network instance (e.g., from PETAce-Network)
     */
    void refill_buffer_sender(const std::shared_ptr<network::Network>& net);

    /**
     * @brief Refill the OT buffer for the receiver
     *
     * @param[in] net The network instance (e.g., from PETAce-Network)
     */
    void refill_buffer_receiver(const std::shared_ptr<network::Network>& net);

    /**
     * @brief Random OT API for the sender
     *
     * @param[in] net The network instance (e.g., from PETAce-Network)
     * @param[in] msgs Messages to send
     * @throws std::invalid_argument if length of LabelType > 128 bits.
     */
    template <typename LabelType>
    void get_random_ot(const std::shared_ptr<network::Network>& net, std::vector<LabelType>& msgs) {
        if (sizeof(LabelType) > 16) {
            throw std::invalid_argument("not support for message which length > 128");
            // todo support message length by using prng
        }
        if (cur_idx_sender_ == buf_sizes_) {
            refill_buffer_sender(net);
            cur_idx_sender_ = 0;
        }
        msgs.resize(2);
        msgs[0] = reinterpret_cast<LabelType*>(&(send_msgs_[cur_idx_sender_][0]))[0];
        msgs[1] = reinterpret_cast<LabelType*>(&(send_msgs_[cur_idx_sender_][1]))[0];

        cur_idx_sender_++;

        return;
    }

    /**
     * @brief Random OT API for the receiver
     *
     * @param[in] net The network instance (e.g., from PETAce-Network)
     * @param[out] choices Choices bits for OT
     * @param[out] msgs Messages received according to the choices
     */
    template <typename LabelType>
    void get_random_ot(const std::shared_ptr<network::Network>& net, OTChoice& choices, LabelType& msgs) {
        if (cur_idx_recver_ == buf_sizes_) {
            refill_buffer_receiver(net);
            cur_idx_recver_ = 0;
        }
        msgs = reinterpret_cast<LabelType*>(&(recv_msgs_[cur_idx_recver_]))[0];
        choices = static_cast<std::int8_t>(verse::bit_from_blocks(ex_choices_, cur_idx_recver_));
        cur_idx_recver_++;
        return;
    }

    /**
     * @brief 2-choose-1 OT API for the sender
     *
     * @param[in] net The network instance (e.g., from PETAce-Network)
     * @param[in] msgs Messages to send
     */
    template <typename MessageType>
    void get_standard_ot(const std::shared_ptr<network::Network>& net, std::vector<std::vector<MessageType>>& msgs) {
        bool* xor_choices = new bool[msgs.size()];
        std::vector<std::array<MessageType, 2>> xor_msgs(msgs.size());
        std::vector<std::vector<MessageType>> ot_msgs(msgs.size());
        for (std::size_t i = 0; i < msgs.size(); i++) {
            get_random_ot(net, ot_msgs[i]);
        }
        recv_bool(net, xor_choices, msgs.size());
        for (std::size_t i = 0; i < msgs.size(); i++) {
            xor_msgs[i][0] = ot_msgs[i][xor_choices[i]] ^ msgs[i][0];
            xor_msgs[i][1] = ot_msgs[i][xor_choices[i] ^ 1] ^ msgs[i][1];
        }
        delete[] xor_choices;

        net->send_data(&xor_msgs[0][0], xor_msgs.size() * sizeof(std::array<MessageType, 2>));

        return;
    }

    /**
     * @brief 2-choose-1 OT API for the receiver
     *
     * @param[in] net The network instance (e.g., from PETAce-Network)
     * @param[in] choices Choices bits for OT
     * @param[out] msgs Messages received
     */
    template <typename MessageType>
    void get_standard_ot(const std::shared_ptr<network::Network>& net, const std::vector<OTChoice>& choices,
            std::vector<MessageType>& msgs) {
        std::vector<OTChoice> rand_choices(choices.size());
        std::vector<MessageType> rand_msgs(choices.size());
        for (std::size_t i = 0; i < choices.size(); i++) {
            get_random_ot(net, rand_choices[i], rand_msgs[i]);
        }

        bool* xor_choices = new bool[choices.size()];
        for (std::size_t i = 0; i < choices.size(); i++) {
            xor_choices[i] = static_cast<bool>(rand_choices[i] ^ choices[i]);
        }
        send_bool(net, xor_choices, choices.size());
        delete[] xor_choices;

        std::vector<std::array<MessageType, 2>> xor_msgs(choices.size());
        net->recv_data(&xor_msgs[0][0], xor_msgs.size() * sizeof(std::array<MessageType, 2>));

        for (std::size_t i = 0; i < choices.size(); i++) {
            msgs.emplace_back(xor_msgs[i][choices[i]] ^ rand_msgs[i]);
        }
        return;
    }

    /**
     * @brief Random Correlated OT API for the sender.
     *
     * @param[in] net The network instance (e.g., from PETAce-Network).
     * @param[in] ot_size The size of correlated ot.
     * @param[in] delta The delta value of correlated ot.
     * @param[in] msgs Messages to send.
     * @throws std::invalid_argument if length of LabelType > 128 bits.
     */
    void get_correlated_ot(const std::shared_ptr<network::Network>& net, const std::size_t ot_size,
            const std::int64_t delta, std::vector<std::vector<std::int64_t>>& msgs);

    /**
     * @brief Random Correlated OT API for the receiver.
     *
     * @param[in] net The network instance (e.g., from PETAce-Network).
     * @param[in] ot_size The size of correlated ot.
     * @param[out] choices Choices bits for OT.
     * @param[out] msgs Messages received according to the choices.
     */
    void get_correlated_ot(const std::shared_ptr<network::Network>& net, const std::size_t ot_size,
            std::vector<std::int8_t>& choices, std::vector<std::int64_t>& msgs);

    /**
     * @brief Random Correlated OT API for the sender.
     *
     * @param[in] net The network instance (e.g., from PETAce-Network).
     * @param[in] ot_size The size of correlated ot.
     * @param[in] delta The delta value of correlated ot.
     * @param[in] msgs Messages to send.
     * @throws std::invalid_argument if length of LabelType > 128 bits.
     */
    void get_batch_cot(const std::shared_ptr<network::Network>& net, const std::size_t ot_size,
            const std::vector<std::int64_t>& delta, std::vector<std::vector<std::int64_t>>& msgs);

    /**
     * @brief Random Correlated OT API for the receiver.
     *
     * @param[in] net The network instance (e.g., from PETAce-Network).
     * @param[in] ot_size The size of correlated ot.
     * @param[out] choices Choices bits for OT.
     * @param[out] msgs Messages received according to the choices.
     */
    void get_batch_cot(const std::shared_ptr<network::Network>& net, const std::size_t ot_size,
            std::vector<std::int8_t>& choices, std::vector<std::int64_t>& msgs);

    /**
     * @brief n-choose-1 OT API for the sender.
     *
     * Refer: Oblivious Transfer and Polynomial Evaluation
     *
     * @param[in] net The network instance (e.g., from PETAce-Network)
     * @param[in] select_size The parameter n
     * @param[in] msgs Messages to send
     */
    void get_nch1_ot(const std::shared_ptr<network::Network>& net, const std::size_t select_size,
            const std::vector<std::vector<std::int64_t>>& msgs);

    /**
     * @brief n-choose-1 OT API for the receiver
     *
     * Refer: Oblivious Transfer and Polynomial Evaluation
     *
     * @param[in] net The network instance (e.g., from PETAce-Network)
     * @param[in] select_size The parameter n
     * @param[in] choices Choices bits for OT
     * @param[out] msgs Messages received
     */
    void get_nch1_ot(const std::shared_ptr<network::Network>& net, const std::size_t select_size,
            const std::vector<std::int8_t>& choices, std::vector<std::int64_t>& msgs);

private:
    std::size_t party_id_ = 0;
    std::size_t buf_sizes_ = 0;
    std::size_t cur_idx_sender_ = 0;
    std::size_t cur_idx_recver_ = 0;
    std::vector<std::array<block, 2>> send_msgs_{};
    std::vector<block> recv_msgs_{};
    std::vector<block> base_choices_{};
    std::vector<block> ex_choices_{};
    std::shared_ptr<solo::PRNG> prng_ = nullptr;
    std::shared_ptr<verse::BaseOtSender> base_ot_sender_ = nullptr;
    std::shared_ptr<verse::BaseOtReceiver> base_ot_recver_ = nullptr;
    std::shared_ptr<verse::OtExtSender> ex_ot_sender_ = nullptr;
    std::shared_ptr<verse::OtExtReceiver> ex_ot_recver_ = nullptr;
};

}  // namespace duet
}  // namespace petace

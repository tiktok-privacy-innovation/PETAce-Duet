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

#include <random>

#include "Eigen/Dense"

#include "network/network.h"

namespace petace {
namespace duet {

inline void send_block(const std::shared_ptr<network::Network>& net, const block* data, std::size_t nblock) {
    net->send_data(data, nblock * sizeof(block));
}

inline void recv_block(const std::shared_ptr<network::Network>& net, block* data, std::size_t nblock) {
    net->recv_data(data, nblock * sizeof(block));
}

inline void send_bool(const std::shared_ptr<network::Network>& net, bool* data, std::size_t length) {
    std::size_t byte_len = (length + 7) / 8;
    std::vector<uint8_t> byte_buf(byte_len, 0);
    for (std::size_t i = 0; i < length; i++) {
        byte_buf[i / 8] |= static_cast<uint8_t>(static_cast<uint8_t>(data[i]) << (7 - (i % 8)));
    }
    net->send_data(byte_buf.data(), byte_len);
}

inline void recv_bool(const std::shared_ptr<network::Network>& net, bool* data, std::size_t length) {
    std::size_t byte_len = (length + 7) / 8;
    std::vector<uint8_t> byte_buf(byte_len);
    net->recv_data(byte_buf.data(), byte_len);
    for (std::size_t i = 0; i < length; i++) {
        data[i] = static_cast<bool>((byte_buf[i / 8] >> (7 - (i % 8))) & 0x1);
    }
}

inline void send_matrix(
        const std::shared_ptr<network::Network>& net, const PlainMatrix<std::int64_t>* data, std::size_t nmatrix) {
    net->send_data(data->data(), nmatrix * data->size() * sizeof(std::int64_t));
}

inline void recv_matrix(
        const std::shared_ptr<network::Network>& net, PlainMatrix<std::int64_t>* data, std::size_t nmatrix) {
    std::vector<std::int64_t> buf(data->size() * nmatrix);
    net->recv_data(&buf[0], nmatrix * data->size() * sizeof(std::int64_t));
    *data = Eigen::Map<PlainMatrix<std::int64_t>>(buf.data(), data->rows(), data->cols());
}

}  // namespace duet
}  // namespace petace
